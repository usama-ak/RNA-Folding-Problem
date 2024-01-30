from utils.file_io import read_pdb, generate_files, generate_csv
from utils.distance_calculation import calculate_distances, distribution_distances, compute_observed_frequency, compute_reference_frequency, compute_pseudo_energy
import os

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_folder = os.path.join(current_dir, '..', 'data')
    train_path = os.path.join(data_folder, 'train-data')
    
    all_distances = []  # Common variable to accumulate distances for all files
    
    for filename in os.listdir(train_path):
        if filename.endswith(".pdb"):
            file_path = os.path.join(train_path, filename)
            
            # Process each file and accumulate distances
            pdb_values = read_pdb(file_path)
            distances = calculate_distances(pdb_values)
            all_distances.extend(distances)

    unique_pairs = set()
    for item in all_distances:
        pair = tuple(sorted([item[0], item[1]]))
        unique_pairs.add(pair)

    distribution = distribution_distances(unique_pairs, all_distances)

    observed_frequency = compute_observed_frequency(distribution)
    reference_frequency = compute_reference_frequency(all_distances)
    scores = compute_pseudo_energy(observed_frequency, reference_frequency)

    generate_files(scores)

    # Specify the CSV filename
    csv_filename = "pseudo_energy_data.csv"
    generate_csv(scores, csv_filename)

    print("Processing complete.")
