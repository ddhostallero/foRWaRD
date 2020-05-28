import numpy as np

# from ProGENI (https://github.com/KnowEnG/ProGENI)
def rwr_vec(node_names, network_matrix, restart_vec, restart_prob, max_iter, tolerance):
    """
    Performs a RWR (Random Walk with Restart) with the given parameters on a
    vector input.
    Input:
        node_names: Name of the nodes in the network
        network_matrix: The probability transition matrix of the network (symmetric)
        restart_vec: The vector representing the restart set
        restart_prob: Probability of restart
        max_iter: Maximum number of iterations for convergence
        tolerance: The threshold used with the residual to determine convergence
    Output:
        num_iter_tmp: Actual number of iterations performed
        residual: The final value of residual
        steady_prob_new: The equlibrium distribution
    """    # Get the number of nodes
    num_nodes = len(node_names)
    no_restart_prob = 1 - restart_prob
    # Compute the initial probability for the nodes
    init_prob = 1/num_nodes
    # Create the vector of probabilities for the nodes
    steady_prob_old = np.empty(num_nodes)
    steady_prob_old.fill(init_prob)
    # Initialize the loop variables (100 is an arbitrary high value)
    residual = 100
    num_iter_tmp = 0
    while (residual > tolerance) and (num_iter_tmp < max_iter):
        steady_prob_new = np.dot(steady_prob_old, network_matrix)
        steady_prob_new *= no_restart_prob
        steady_prob_new += restart_prob * restart_vec
        # Calculate the residual -- the sum of the absolute
        # differences of the new node probability vector minus the old
        residual = abs(steady_prob_new - steady_prob_old).sum()
        num_iter_tmp += 1
        steady_prob_old = steady_prob_new.copy()
    return (num_iter_tmp, residual, steady_prob_new)
