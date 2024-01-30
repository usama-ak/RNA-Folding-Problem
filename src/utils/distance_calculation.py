import math
from math import sqrt, pow

def calculate_distances(values_list):
    """
    Calculates intrachain distances between C3' atoms.

    Parameters:
    - values_list: List of lists containing information about C3' atoms.

    Returns:
    - List of lists containing intrachain distances.
    """
    distances_list = []

    for i in range(len(values_list)):
        for j in range(i + 4, len(values_list)):
            # Check if residues belong to the same chain before calculating distance
            if values_list[i][2] == values_list[j][2]:
                d = sqrt(pow(float(values_list[i][4]) - float(values_list[j][4]), 2) +
                         pow(float(values_list[i][5]) - float(values_list[j][5]), 2) +
                         pow(float(values_list[i][6]) - float(values_list[j][6]), 2))
                if d <= 20:
                    distances_list.append([values_list[i][1], values_list[j][1], d])

    return distances_list

def distribution_distances(unique_pairs, distances):
    """
    Groups distances by unique pairs.

    Parameters:
    - unique_pairs: Set of unique residue pairs.
    - distances: List of lists containing intrachain distances.

    Returns:
    - Dictionary associating each unique pair with its distances.
    """
    associations = {}
    for pair in unique_pairs:
        associated_distances = []
        for item in distances:
            if set(pair) == set([item[0], item[1]]):
                associated_distances.append(item[2])
        associations[pair] = associated_distances
    return associations

def compute_observed_frequency(associations):
    """
    Computes observed frequency of distances within intervals.

    Parameters:
    - associations: Dictionary associating unique pairs with their distances.

    Returns:
    - Dictionary associating unique pairs with observed frequency intervals.
    """
    observed_frequency = {}
    for pair, distances in associations.items():
        # Initialize a list with 20 intervals
        interval_counts = [0] * 20
        # Count the occurrences of each distance within the intervals
        for distance in distances:
            interval_index = int(distance / 1)
            interval_counts[interval_index] += 1
        # Compute the observed frequency for each interval
        observed_frequency[pair] = [count / len(distances) if len(distances) > 0 else 0 for count in interval_counts]
    
    return observed_frequency

def compute_reference_frequency(distances_list):
    """
    Computes reference frequency of all distances within intervals.

    Parameters:
    - distances_list: List of lists containing all distances.

    Returns:
    - List of reference frequency intervals.
    """
    all_distances = [item[2] for item in distances_list]
    interval_counts = [0] * 20
    for distance in all_distances:
        interval_index = int(distance / 1)
        interval_counts[interval_index] += 1
    reference_frequency = [count / len(all_distances) if len(all_distances) > 0 else 0 for count in interval_counts]
    return reference_frequency


def compute_pseudo_energy(observed_frequency, reference_frequency):
    """
    Computes pseudo-energy scores based on observed and reference frequencies.

    Parameters:
    - observed_frequency: Dictionary associating unique pairs with observed frequency intervals.
    - reference_frequency: List of reference frequency intervals.

    Returns:
    - Dictionary associating unique pairs with pseudo-energy scores.
    """
    pseudo_energy_scores = {}

    for pair, observed_interval in observed_frequency.items():
        reference_interval = reference_frequency  

        # Compute pseudo-energy scores for each interval
        scores = []
        for obs, ref in zip(observed_interval, reference_interval):
            if ref != 0 and obs != 0:
                scores.append(-math.log(obs / ref))
            else:
                scores.append(float('inf'))
        
        pseudo_energy_scores[pair] = scores

    return pseudo_energy_scores
