from importlib import resources
import itertools
from tqdm import tqdm
import numpy as np

from . import data


class Msa:
    def load_blosum(self) -> dict[str, dict[str, int]]:
        lines = resources.read_text(data, "blosum62.txt").splitlines()

        blosum: dict[str, dict[str, int]] = {}
        amino_acids = []
        for line in lines:
            if line[0] == " ":
                amino_acids = line.split()
                for amino_acid_j in amino_acids:
                    blosum[amino_acid_j] = {}
            else:
                amino_acid_i = line[0]
                values = line[1:].split()
                for i, amino_acid_j in enumerate(amino_acids):
                    blosum[amino_acid_i][amino_acid_j] = int(values[i])

        return blosum

    def align(self, sequences: list[str]):
        blosum = self.load_blosum()
        lengths = [len(sequence) for sequence in sequences]
        dim = tuple([l + 1 for l in lengths])
        n = len(sequences)
        dp = np.full(dim, -10, dtype=int)

        # fill dp
        for indices in tqdm(np.ndindex(dim), total=int(np.prod(dim))):
            for delta in itertools.product([0, 1], repeat=n):
                if delta == (0,) * n:
                    continue
                if any([indices[i] - delta[i] < 0 for i in range(n)]):
                    continue

                score = 0
                indices_to_use = np.where(np.array(delta) == 1)[0]
                for i, j in itertools.combinations(np.arange(n), 2):
                    if indices[i] >= lengths[i] or indices[j] >= lengths[j]:
                        continue
                    amino_acid_i = "*"
                    amino_acid_j = "*"
                    if i in indices_to_use:
                        amino_acid_i = sequences[i][indices[i]]
                    if j in indices_to_use:
                        amino_acid_j = sequences[j][indices[j]]
                    score += blosum[amino_acid_i][amino_acid_j]

                indices_delta = tuple(np.array(indices) - np.array(delta))
                dp[indices] = max(dp[indices], score + dp[indices_delta])

        # traceback
        alignments: list[list[str]] = []
        indices = tuple(lengths)
        while np.any(indices):
            max_score = -np.inf
            max_delta = (-1,) * n
            for delta in itertools.product([0, 1], repeat=n):
                if delta == (0,) * n:
                    continue
                if any([indices[i] - delta[i] < 0 for i in range(n)]):
                    continue

                indices_delta = tuple(np.array(indices) - np.array(delta))
                max_score = max(max_score, dp[indices_delta])
                if dp[indices_delta] == max_score:
                    max_delta = delta

            indices = tuple(np.array(indices) - np.array(max_delta))
            alignments = [
                [
                    sequences[i][indices[i]] if max_delta[i] == 1 else "-"
                    for i in range(n)
                ]
            ] + alignments

        result = []
        for i in range(n):
            result.append("".join([alignment[i] for alignment in alignments]))
        return result
