import math

def read_scores(score_file):
    """
    Reads scoring values from a file.

    Parameters:
    - score_file: Path to the file containing the scores.

    Returns:
    - List of scores for intervals of distances between C3' atoms.
    """
    scores = []
    with open(score_file, "r") as file:
        for line in file:
            score = float(line.strip())
            scores.append(score)
    return scores

def linear_interpolation(distance, scores):
    """
    Perform linear interpolation to calculate a score based on the distance.

    Parameters:
    - distance: Distance between C3' atoms.
    - scores: List of scores for intervals of distances between C3' atoms.

    Returns:
    - Score computed using linear interpolation.
    """
    interval_index = int(distance)
    if interval_index < 0 or interval_index >= len(scores):
        return float('inf')  # Distance falls outside the valid range

    # Handle the upper bound of the range separately
    if interval_index == len(scores) - 1:
        return scores[-1]

    # Calculate the fractional distance within the interval
    fractional_distance = distance - interval_index

    # Compute linear interpolation
    score_low = scores[interval_index]
    score_high = scores[interval_index + 1]
    interpolated_score = score_low + fractional_distance * (score_high - score_low)
    return interpolated_score

def calculate_gibbs_free_energy(scores):
    """
    Calculates the estimated Gibbs free energy based on scores.

    Parameters:
    - scores: List of scores for distances between C3' atoms.

    Returns:
    - The estimated Gibbs free energy.
    """
    # Filter out 'NaN' values
    scores = [score for score in scores if not math.isnan(score)]
    
    # Calculate sum of non-'NaN' scores
    return sum(scores)
