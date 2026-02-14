"""Ensemble analysis: RMSD computation, clustering, and consensus detection."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

logger = logging.getLogger("rnapipey")

# Backbone atom names used for structural comparison
BACKBONE_ATOMS = {"C3'", "P"}


@dataclass
class ClusterInfo:
    cluster_id: int
    members: list[Path]
    member_predictors: list[str]
    representative: Path
    mean_rmsd: float
    is_consensus: bool  # members from 2+ different predictors


@dataclass
class EnsembleResult:
    clusters: list[ClusterInfo]
    rmsd_matrix: np.ndarray
    pdb_files: list[Path]
    predictor_labels: list[str]


def _extract_backbone_coords(structure_path: Path) -> np.ndarray:
    """Extract C3' and P atom coordinates from a PDB or CIF file."""
    suffix = structure_path.suffix.lower()

    if suffix == ".cif":
        from Bio.PDB import MMCIFParser
        parser = MMCIFParser(QUIET=True)
    else:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("s", str(structure_path))
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() in BACKBONE_ATOMS:
                        coords.append(atom.get_vector().get_array())
        break  # only first model

    if not coords:
        raise ValueError(f"No backbone atoms (C3'/P) found in {structure_path}")
    return np.array(coords)


def compute_pairwise_rmsd(pdb_files: list[Path]) -> np.ndarray:
    """Compute pairwise RMSD matrix over backbone atoms (C3', P).

    Uses BioPython Superimposer for optimal superposition.
    """
    from Bio.PDB import Superimposer

    n = len(pdb_files)
    rmsd_matrix = np.zeros((n, n))

    # Pre-extract coordinates
    all_coords: list[np.ndarray] = []
    for f in pdb_files:
        try:
            coords = _extract_backbone_coords(f)
            all_coords.append(coords)
        except Exception as e:
            logger.warning("Could not parse %s: %s", f, e)
            all_coords.append(np.array([]))

    sup = Superimposer()
    for i in range(n):
        for j in range(i + 1, n):
            ci, cj = all_coords[i], all_coords[j]
            if ci.size == 0 or cj.size == 0:
                rmsd_matrix[i, j] = rmsd_matrix[j, i] = float("inf")
                continue

            # Truncate to common length (different predictors may output
            # slightly different numbers of atoms)
            min_len = min(len(ci), len(cj))
            if min_len < 3:
                rmsd_matrix[i, j] = rmsd_matrix[j, i] = float("inf")
                continue

            try:
                sup.set_atoms(
                    _coords_to_atoms(ci[:min_len]),
                    _coords_to_atoms(cj[:min_len]),
                )
                rmsd_val = sup.rms
                rmsd_matrix[i, j] = rmsd_matrix[j, i] = rmsd_val
            except Exception:
                rmsd_matrix[i, j] = rmsd_matrix[j, i] = float("inf")

    return rmsd_matrix


def _coords_to_atoms(coords: np.ndarray) -> list:
    """Convert coordinate array to pseudo-Atom objects for Superimposer."""
    from Bio.PDB.Atom import Atom
    from Bio.PDB.vectors import Vector

    atoms = []
    for i, c in enumerate(coords):
        a = Atom(
            name="CA",
            coord=c,
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=" CA ",
            serial_number=i,
            element="C",
        )
        atoms.append(a)
    return atoms


def cluster_structures(
    pdb_files: list[Path],
    predictor_labels: list[str],
    cutoff: float = 5.0,
) -> EnsembleResult:
    """Cluster structures by RMSD using hierarchical clustering.

    Args:
        pdb_files: List of structure file paths (PDB or CIF).
        predictor_labels: Predictor name for each file (same length as pdb_files).
        cutoff: RMSD cutoff in Angstroms for cluster formation.

    Returns:
        EnsembleResult with cluster assignments, RMSD matrix, etc.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    n = len(pdb_files)
    if n <= 1:
        # Single structure: trivial cluster
        clusters = []
        if n == 1:
            clusters.append(ClusterInfo(
                cluster_id=1,
                members=[pdb_files[0]],
                member_predictors=[predictor_labels[0]],
                representative=pdb_files[0],
                mean_rmsd=0.0,
                is_consensus=False,
            ))
        return EnsembleResult(
            clusters=clusters,
            rmsd_matrix=np.zeros((n, n)),
            pdb_files=pdb_files,
            predictor_labels=predictor_labels,
        )

    logger.info("Computing pairwise RMSD for %d structures...", n)
    rmsd_matrix = compute_pairwise_rmsd(pdb_files)

    # Replace inf with large value for clustering
    finite_max = np.max(rmsd_matrix[np.isfinite(rmsd_matrix)]) if np.any(np.isfinite(rmsd_matrix)) else 100.0
    rmsd_safe = np.where(np.isfinite(rmsd_matrix), rmsd_matrix, finite_max * 2)

    # Hierarchical clustering (average linkage)
    condensed = squareform(rmsd_safe)
    Z = linkage(condensed, method="average")
    labels = fcluster(Z, t=cutoff, criterion="distance")

    # Build ClusterInfo for each cluster
    cluster_ids = sorted(set(labels))
    clusters: list[ClusterInfo] = []

    for cid in cluster_ids:
        member_indices = [i for i, l in enumerate(labels) if l == cid]
        members = [pdb_files[i] for i in member_indices]
        preds = [predictor_labels[i] for i in member_indices]

        # Mean intra-cluster RMSD
        if len(member_indices) > 1:
            sub_rmsds = []
            for ii, a in enumerate(member_indices):
                for b in member_indices[ii + 1:]:
                    sub_rmsds.append(rmsd_matrix[a, b])
            mean_rmsd = float(np.mean(sub_rmsds))
        else:
            mean_rmsd = 0.0

        # Representative: closest to centroid (lowest average RMSD to other members)
        if len(member_indices) > 1:
            avg_rmsds = []
            for i in member_indices:
                avg = np.mean([rmsd_matrix[i, j] for j in member_indices if j != i])
                avg_rmsds.append(avg)
            rep_idx = member_indices[int(np.argmin(avg_rmsds))]
        else:
            rep_idx = member_indices[0]

        is_consensus = len(set(preds)) >= 2

        clusters.append(ClusterInfo(
            cluster_id=int(cid),
            members=members,
            member_predictors=preds,
            representative=pdb_files[rep_idx],
            mean_rmsd=mean_rmsd,
            is_consensus=is_consensus,
        ))

    # Sort by cluster size descending
    clusters.sort(key=lambda c: len(c.members), reverse=True)

    logger.info(
        "Clustering complete: %d clusters from %d structures (cutoff=%.1f A)",
        len(clusters), n, cutoff,
    )
    consensus_count = sum(1 for c in clusters if c.is_consensus)
    if consensus_count:
        logger.info("  %d consensus cluster(s) (2+ predictors agree)", consensus_count)

    return EnsembleResult(
        clusters=clusters,
        rmsd_matrix=rmsd_matrix,
        pdb_files=pdb_files,
        predictor_labels=predictor_labels,
    )
