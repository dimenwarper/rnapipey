"""Batch inference script for RhoFold+: loads model once, runs multiple seeds.

Invoked as a subprocess (like inference.py) but processes all seeds in a
single process to avoid repeated model loading (~9 min per load).

Usage:
    python rhofold_batch.py \
        --input_fas query.fasta \
        --ckpt /path/to/rhofold_pretrained_params.pt \
        --seeds 0,1,2,3,4 \
        --output_base_dir /path/to/output \
        --device cuda:0 \
        --single_seq_pred True
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from huggingface_hub import snapshot_download

from rhofold.config import rhofold_config
from rhofold.rhofold import RhoFold
from rhofold.utils import get_device, save_ss2ct, timing
from rhofold.utils.alphabet import get_features


logger = logging.getLogger("RhoFold-Batch")


def setup_logger(output_dir: str) -> None:
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")

    fh = logging.FileHandler(f"{output_dir}/batch_log.txt", mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.DEBUG)
    sh.setFormatter(formatter)
    logger.addHandler(sh)


@torch.no_grad()
def main() -> None:
    parser = argparse.ArgumentParser(description="RhoFold+ batch inference")
    parser.add_argument("--input_fas", required=True)
    parser.add_argument("--ckpt", default="./pretrained/rhofold_pretrained_params.pt")
    parser.add_argument("--device", default=None)
    parser.add_argument("--single_seq_pred", default="False")
    parser.add_argument("--input_a3m", default=None)
    parser.add_argument("--seeds", required=True, help="Comma-separated seeds, e.g. 0,1,2,3,4")
    parser.add_argument("--output_base_dir", required=True)
    parser.add_argument("--mc_dropout", action="store_true",
                        help="Re-enable dropout at inference for MC Dropout diversity")
    parser.add_argument("--noise_scale", type=float, default=0.0,
                        help="Gaussian noise scale to add to post-embedding features (0 = off)")

    args = parser.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    output_base = Path(args.output_base_dir)
    output_base.mkdir(parents=True, exist_ok=True)

    setup_logger(str(output_base))

    # ---- Load model ONCE ----
    logger.info("Constructing RhoFold")
    model = RhoFold(rhofold_config)

    ckpt_path = Path(args.ckpt)
    if not ckpt_path.exists():
        logger.info("Downloading checkpoint to %s", ckpt_path.parent)
        try:
            snapshot_download(repo_id="cuhkaih/rhofold", local_dir=ckpt_path.parent)
        except Exception as e:
            logger.error("Could not download checkpoint: %s", e)
            raise
    else:
        logger.info("Checkpoint already exists at %s, skipping download", ckpt_path)

    logger.info("Loading checkpoint %s", args.ckpt)
    model.load_state_dict(
        torch.load(args.ckpt, map_location=torch.device("cpu"))["model"]
    )
    model.eval()

    # ---- MC Dropout: re-enable dropout layers in eval mode ----
    if args.mc_dropout:
        dropout_count = 0
        for module in model.modules():
            if isinstance(module, nn.Dropout):
                module.train()
                dropout_count += 1
        logger.info("MC Dropout enabled: set %d Dropout layers to train mode", dropout_count)

    # ---- Move to device ONCE ----
    device = get_device(args.device)
    logger.info("Using device %s", device)
    model = model.to(device)

    # ---- Prepare input features ONCE ----
    single_seq = str(args.single_seq_pred).lower() in ("true", "1", "yes")
    input_a3m = args.input_a3m
    if single_seq:
        input_a3m = args.input_fas
        logger.info("Single-sequence mode: using input_fas as MSA")
    elif input_a3m:
        logger.info("Using provided MSA: %s", input_a3m)
    else:
        logger.info("No MSA provided and single_seq_pred=False; using input_fas as fallback")
        input_a3m = args.input_fas

    data_dict = get_features(args.input_fas, input_a3m)

    # ---- Embedding noise hook ----
    # Noise must be added AFTER nn.Embedding (which requires integer tokens),
    # so we use a forward hook on the MSA embedder to perturb its output.
    noise_scale = args.noise_scale
    _noise_active = False  # mutable flag toggled per seed

    def _embedding_noise_hook(module, input, output):
        """Add Gaussian noise to (msa_fea, pair_fea) after embedding."""
        if not _noise_active:
            return output
        msa_fea, pair_fea = output
        msa_fea = msa_fea + torch.randn_like(msa_fea) * noise_scale
        pair_fea = pair_fea + torch.randn_like(pair_fea) * noise_scale
        return msa_fea, pair_fea

    hook_handle = None
    if noise_scale > 0:
        hook_handle = model.msa_embedder.register_forward_hook(_embedding_noise_hook)
        logger.info("Registered embedding noise hook (scale=%.4f)", noise_scale)

    # ---- Loop over seeds ----
    n_seeds = len(seeds)
    use_mc_dropout = args.mc_dropout

    for i, seed in enumerate(seeds):
        os.environ["PYTHONHASHSEED"] = str(seed)
        is_first_seed = (i == 0)

        run_dir = output_base / f"run_{seed}"
        run_dir.mkdir(parents=True, exist_ok=True)

        # Seed 0 (first seed) is always vanilla: no dropout, no noise
        if is_first_seed and (use_mc_dropout or noise_scale > 0):
            logger.info(
                "Seed %d (%d/%d) — VANILLA baseline (no dropout/noise) — output: %s",
                seed, i + 1, n_seeds, run_dir,
            )
            # Temporarily disable dropout for vanilla run
            if use_mc_dropout:
                for module in model.modules():
                    if isinstance(module, nn.Dropout):
                        module.eval()
            _noise_active = False
        else:
            logger.info("Seed %d (%d/%d) — output: %s", seed, i + 1, n_seeds, run_dir)
            if not is_first_seed and noise_scale > 0:
                torch.manual_seed(seed)
                _noise_active = True
                logger.info("Embedding noise active (scale=%.4f, seed=%d)", noise_scale, seed)
            else:
                _noise_active = False

        with timing(f"Seed {seed} inference", logger=logger):
            outputs = model(
                tokens=data_dict["tokens"].to(device),
                rna_fm_tokens=data_dict["rna_fm_tokens"].to(device),
                seq=data_dict["seq"],
            )

        # Re-enable dropout after vanilla run
        if is_first_seed and use_mc_dropout:
            for module in model.modules():
                if isinstance(module, nn.Dropout):
                    module.train()

        output = outputs[-1]

        # Secondary structure (.ct)
        ss_prob_map = torch.sigmoid(output["ss"][0, 0]).data.cpu().numpy()
        save_ss2ct(ss_prob_map, data_dict["seq"], f"{run_dir}/ss.ct", threshold=0.5)

        # Distance maps + pLDDT (.npz)
        np.savez_compressed(
            f"{run_dir}/results.npz",
            dist_n=torch.softmax(output["n"].squeeze(0), dim=0).data.cpu().numpy(),
            dist_p=torch.softmax(output["p"].squeeze(0), dim=0).data.cpu().numpy(),
            dist_c=torch.softmax(output["c4_"].squeeze(0), dim=0).data.cpu().numpy(),
            ss_prob_map=ss_prob_map,
            plddt=output["plddt"][0].data.cpu().numpy(),
        )

        # 3D structure (.pdb)
        node_cords_pred = output["cord_tns_pred"][-1].squeeze(0)
        model.structure_module.converter.export_pdb_file(
            data_dict["seq"],
            node_cords_pred.data.cpu().numpy(),
            path=f"{run_dir}/unrelaxed_model.pdb",
            chain_id=None,
            confidence=output["plddt"][0].data.cpu().numpy(),
            logger=logger,
        )

        logger.info("Seed %d (%d/%d) done", seed, i + 1, n_seeds)

    if hook_handle is not None:
        hook_handle.remove()

    logger.info("Batch inference complete: %d seeds processed", n_seeds)


if __name__ == "__main__":
    main()
