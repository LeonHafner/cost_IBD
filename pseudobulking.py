#!/usr/bin/env python3
import argparse
from typing import Callable, Sequence, Union
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from functools import reduce
from operator import and_


def pseudobulk(
    adata: AnnData,
    *,
    groupby: Union[str, Sequence[str]],
    aggr_fun: Callable = np.sum,
    min_obs: int = 0,
) -> AnnData:
    """
    Pseudobulk by summing counts per dataset × sample × cell type.
    Keeps integer counts for DESeq2.
    """
    if isinstance(groupby, str):
        groupby = [groupby]

    missing = [g for g in groupby if g not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing required obs columns: {missing}")

    # unique combinations to iterate
    combos = adata.obs.loc[:, groupby].drop_duplicates()

    # precompute boolean masks as NumPy arrays
    masks = {
        col: {
            val: (adata.obs[col].values == val)  # NumPy bool array
            for val in combos[col].unique()
        }
        for col in groupby
    }

    expr_agg = []
    obs_rows = []

    for comb in combos.itertuples(index=False):
        mask = reduce(and_, (masks[col][val] for col, val in zip(groupby, comb)))
        n = int(mask.sum())
        if n < min_obs:
            continue

        X_sub = adata.X[mask, :]

        # aggregate
        try:
            row = aggr_fun(X_sub, axis=0).A1  # sparse path
        except AttributeError:
            row = aggr_fun(X_sub, axis=0)

        # ensure plain NumPy array
        row = np.asarray(row)

        # if aggregation is sum, enforce integer dtype for DESeq2
        if aggr_fun is np.sum:
            # safe cast after rounding to protect against float sparse sums
            row = np.rint(row).astype(np.int64, copy=False)

        obs_row = dict(zip(groupby, comb))
        obs_row["n_obs"] = n
        expr_agg.append(row)
        obs_rows.append(obs_row)

    if not expr_agg:
        raise ValueError("No groups passed min_obs. Output would be empty.")

    X = np.vstack(expr_agg)

    return AnnData(
        X=X,
        var=adata.var.copy(),
        obs=pd.DataFrame.from_records(obs_rows),
        uns=adata.uns.copy(),
    )


def main():
    parser = argparse.ArgumentParser(
        prog="Pseudobulking",
        description="Aggregate cells into pseudobulks by dataset, sample, and annotation:coarse",
    )
    parser.add_argument("--input", required=True, help="Input .h5ad")
    parser.add_argument("--output", required=True, help="Output .h5ad")
    parser.add_argument("--min-cells", type=int, default=50, help="Min cells per group")
    parser.add_argument(
        "--groupby",
        nargs="*",
        default=["dataset", "sample", "annotation:coarse"],
        help="obs columns to group by",
    )
    args = parser.parse_args()

    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Groupby: {args.groupby}")
    print(f"Min cells: {args.min_cells}")

    adata = sc.read_h5ad(args.input)

    adata_pb = pseudobulk(
        adata=adata,
        groupby=args.groupby,
        aggr_fun=np.sum,      # sums give counts for DESeq2
        min_obs=args.min_cells,
    )

    # add a compact group label
    adata_pb.obs["pseudobulk_id"] = (
        adata_pb.obs[args.groupby].astype(str).agg("|".join, axis=1)
    )

    adata_pb.write_h5ad(args.output)
    print("Done.")


if __name__ == "__main__":
    main()
