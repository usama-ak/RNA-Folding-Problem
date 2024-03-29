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

def GetScore(distance, scores):
    """
    Perform linear interpolation to calculate a score based on the distance.

    Parameters:
    - distance: Distance between C3' atoms.
    - scores: List of scores for intervals of distances between C3' atoms.

    Returns:
    - Score corresponding to the interval where the distance falls.
    """
    interval_index = int(distance)
    if interval_index < 0 or interval_index >= len(scores):
        return float('inf')  ## If the distance falls outside the valid range, the obtained value is equal to inf.

    return scores[interval_index]


def calculate_gibbs_free_energy(scores):
    """
    Calculates the estimated Gibbs free energy based on scores.

    Parameters:
    - scores: List of scores for distances between C3' atoms.

    Returns:
    - The estimated Gibbs free energy.
    """
    ## We first of all filter out 'NaN' and 'inf' values 
    scores = [score for score in scores if not (math.isnan(score) or math.isinf(score))]
    
    ## Once we filtered the Nan and inf values, we calculate sum of non-'NaN' and non-'inf' scores
    return sum(scores)
