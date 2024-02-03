import os
import argparse
from src.utils.distance_calculation import calculate_distances
from src.utils.file_io import read_pdb
from src.utils.score_calculation import read_scores, GetScore, calculate_gibbs_free_energy

## This class has been created to estimate GibbsFreeEnergy of PDB files coming from the PDB database
class GibbsFreeEnergyCalculator:
    def __init__(self, pdb_file_path, scores_dir="data/scores/"):
        self.pdb_file_path = pdb_file_path
        self.scores_dir = scores_dir

        ## Here, we map nucleotide pairs to score file names
        self.score_files = {
            ('A', 'A'): "A_A_scores.txt",
            ('A', 'U'): "A_U_scores.txt",
            ('A', 'G'): "A_G_scores.txt",
            ('A', 'C'): "A_C_scores.txt",
            ('U', 'U'): "U_U_scores.txt",
            ('U', 'G'): "G_U_scores.txt",
            ('U', 'C'): "C_U_scores.txt",
            ('G', 'G'): "G_G_scores.txt",
            ('G', 'C'): "C_G_scores.txt",
            ('C', 'C'): "C_C_scores.txt",
        }
        ## Then, we append reverse pairs
        self.score_files.update({pair[::-1]: filename for pair, filename in self.score_files.items()})


    def calculate_gibbs_free_energy(self):
        ## We then process the distances and initialize a dictionary to store scores for each pair
        pdb_values = read_pdb(self.pdb_file_path)
        distances = calculate_distances(pdb_values)

        all_scores = {}

        ## A for loop is created to iterate over distances and calculate scores
        for distance in distances:
            nucleotide_1, nucleotide_2, distance_value = distance
            distance_key = tuple(distance)

            ## The appropriate score file is obtained based on nucleotide pair
            score_file = self.score_files.get((nucleotide_1, nucleotide_2))

            ## Constructing the full path to the score file
            full_score_file_path = os.path.join(self.scores_dir, score_file)

            ## Then, we read scores from the score file
            scores = read_scores(full_score_file_path)

            ## Calculating score for current distance using linear interpolation and storing the score on "score" variable
            score = GetScore(distance_value, scores)
            
            all_scores[distance_key] = score

        ## Finally, we calculate Gibbs free energy
        gibbs_free_energy = calculate_gibbs_free_energy(all_scores.values())

        return gibbs_free_energy

if __name__ == "__main__":
    ## Parsing arguments are defined in this section
    parser = argparse.ArgumentParser(description="Calculate Gibbs free energy for RNA structures")
    parser.add_argument("--pdb-file", help="Path to the PDB file containing the RNA structure")
    parser.add_argument("--scores-dir", default="data/scores/", help="Directory containing score files")
    args = parser.parse_args()

    if args.pdb_file:
        ## Use the specified PDB file
        calculator = GibbsFreeEnergyCalculator(args.pdb_file, args.scores_dir)
        gibbs_free_energy = calculator.calculate_gibbs_free_energy()
        print(f"Estimated Gibbs Free Energy: {gibbs_free_energy}")
    else:
        ## No PDB file specified, process example files from the examples directory
        examples_dir = "data/examples/"
        for file_name in os.listdir(examples_dir):
            if file_name.endswith(".pdb"):
                example_file_path = os.path.join(examples_dir, file_name)
                calculator = GibbsFreeEnergyCalculator(example_file_path, args.scores_dir)
                gibbs_free_energy = calculator.calculate_gibbs_free_energy()
                print(f"Gibbs Free Energy for {file_name}: {gibbs_free_energy}")


