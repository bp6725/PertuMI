# PertuMI

Active Model Identification(AMI) with FDR guarantees.

AMI algorithms receive a gene-gene network and genes scores as input and report subnetworks that show significant over-representation of high-scored genes. Then, those subnetworks are compared to known pathways.

Unforthently, those algorithms are prone to high rates of false discovery. Studies showed that most algorithms return similar enrichment scores for both original and completely randomized data (permutated data).

This project uses a bandit allocation-based permutation test for AMI with false discovery rate(FDR) guarantees.

