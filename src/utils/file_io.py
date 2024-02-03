import os
import csv

def read_pdb(file_path):
    """
    Reads a PDB file and extracts relevant information for C3' atoms.

    Parameters:
    - file_path: Path to the PDB file.

    Returns:
    - List of lists containing information about C3' atoms.
    """
    values_list = []
    model_found = False
    ## Here, we open a pdb file
    with open(file_path, "r") as pdbfile:
        for line in pdbfile:
            if line.startswith("MODEL") and "1" in line:
                model_found = True
            elif model_found and line.startswith("ATOM"):
                ## Then we define column positions to extract information that we need
                column_positions = [
                    (11, 18),    ## atom
                    (19, 20),    ## nucleotide
                    (21, 22),    ## chain
                    (23, 26),    ## residue number
                    (27, 37),    ## X
                    (38, 45),    ## Y
                    (46, 53)     ## Z
                ]
                values = [line[start:end].strip() for start, end in column_positions]
                if values[0] == "C3'":
                    values_list.append(values)
            elif model_found and line.startswith("ENDMDL"): ## In case the file contains many models
                break  ## We add a break condition to stop reading the pdb file at the end of the first model
            elif not model_found and line.startswith("ATOM"): ## If it is not the case :
                ## We once again define column positions for relevant information
                column_positions = [
                    (11, 18),    ## atom
                    (19, 20),    ## nucleotide
                    (21, 22),    ## chain
                    (23, 26),    ## residue number
                    (27, 37),    ## X
                    (38, 45),    ## Y
                    (46, 53)     ## Z
                ]
                values = [line[start:end].strip() for start, end in column_positions]
                if values[0] == "C3'":
                    values_list.append(values)

    return values_list

def generate_files(pseudo_energy_scores):
    """
    Generates text files containing pseudo-energy scores.

    Parameters:
    - pseudo_energy_scores: Dictionary associating unique pairs with pseudo-energy scores.

    Returns:
    - None
    """
    for pair, scores in pseudo_energy_scores.items():
        file_name = f"{pair[0]}_{pair[1]}_scores.txt"
        with open(file_name, "w") as file:
            for score in scores:
                file.write(f"{score:.4f}\n")
        print(f"File '{file_name}' created.")

def generate_csv(pseudo_energy_scores, csv_filename):
    """
    Generates a CSV file containing distances and pseudo-energy scores.

    Parameters:
    - pseudo_energy_scores: Dictionary associating unique pairs with pseudo-energy scores.
    - csv_filename: Name of the CSV file to be generated.

    Returns:
    - None
    """
    with open(csv_filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        ## Header is written on the csv file
        csv_writer.writerow(['Nucleotide 1', 'Nucleotide 2'] + [f'Score_{i}' for i in range(1, 21)])
        ## Then the data is written on the same file
        for pair, scores in pseudo_energy_scores.items():
            nucleotide_1, nucleotide_2 = pair
            scores = scores[:]
            csv_writer.writerow([nucleotide_1, nucleotide_2] + scores)
    print("CSV file generated.")
