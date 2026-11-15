"""
Tests for utils.from_rbh, which reads linkage-group sizes and the entropy value m from a .rbh file.

Run from the repository root with: python -m pytest tests/
"""

import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'inversion_sim'))

import utils  # noqa: E402

EXAMPLE_RBH = os.path.join(REPO, 'example_scripts', 'HVU_RES_xy_reciprocal_best_hits.coloredby_BCnS_LGs.plotted.rbh')


def write_rbh(path):
    """
    write a small .rbh file: on scaffold ORG1, group G1 has 3 genes and group G2 has 2 genes, in the order A B A B A;
    a G1 gene on another scaffold and a gene of a third group must not be counted
    """

    rows = [
        ('G1', 'ORG1', 100),
        ('G2', 'ORG1', 200),
        ('G1', 'ORG1', 300),
        ('G2', 'ORG1', 400),
        ('G1', 'ORG1', 500),
        ('G1', 'ORG2', 600),
        ('G3', 'ORG1', 700),
    ]
    with open(path, 'w') as f:
        f.write('rbh\tgene_group\tORG_scaf\tORG_pos\tother_column\n')
        for i, (group, scaf, pos) in enumerate(rows):
            f.write('rbh_{}\t{}\t{}\t{}\tx\n'.format(i, group, scaf, pos))


def test_sizes_are_gene_counts(tmp_path):
    rbh = tmp_path / 'small.rbh'
    write_rbh(rbh)

    Asize, Bsize, m = utils.from_rbh(str(rbh), 'G1', 'G2', 'ORG1', 'ORG')

    assert (Asize, Bsize) == (3, 2)
    # ABABA has 4 transitions: m = (4-1)/(2*3*2/(3+2)-1)
    assert abs(m - 3/(12/5-1)) < 1e-9


def test_example_rbh_sizes():
    # 65 B1 genes and 44 B2 genes are on HVU7 in the example .rbh file
    Asize, Bsize, m = utils.from_rbh(EXAMPLE_RBH, 'B1', 'B2', 'HVU7', 'HVU')

    assert (Asize, Bsize) == (65, 44)
    assert round(m, 3) == 1.049
