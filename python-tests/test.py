import pysps
import pydoc
import numpy as np
import random


# A dummy function for testing
def calc_objective_func(x: np.ndarray) -> float:
    if True:
        # Simple distance
        return -np.linalg.norm(x - 0.2)
    else:
        # Modified Rosenbrock
        num_dims = x.shape[0]
        sum = 0.0
        for dim in range(num_dims - 1):
            sum -= 100.0 * (4.0 * x[dim + 1] - 16.0 * x[dim] * x[dim]) * (4.0 * x[dim + 1] - 16.0 * x[dim] * x[dim])
            sum -= (1.0 - 4.0 * x[dim]) * (1.0 - 4.0 * x[dim])
        return sum


if __name__ == "__main__":
    print(pydoc.plain(pydoc.render_doc(pysps)))
    np.set_printoptions(precision=3)
    pysps.set_seed(random.randint(0, 65535))

    num_candidates = 5
    num_iters = 20
    num_dims = 8
    depth = 12
    use_map_hyperparams = True
    radius = int((num_candidates - 1) / 2)

    optimizer = pysps.Optimizer(num_dims, use_map_hyperparams)

    for iter in range(num_iters):

        plane = optimizer.retrieve_search_plane()

        prev_cell_indices = []
        x_best = None

        for d in range(depth):

            f_best = None
            cell_index_best = None

            for i in range(num_candidates):
                for j in range(num_candidates):
                    cell_index = (i - radius, j - radius)

                    x = plane.calc_grid_parameters(cell_index, num_candidates, 0.5, prev_cell_indices)
                    f = calc_objective_func(x)

                    if not f_best or f_best < f:
                        x_best = x
                        f_best = f
                        cell_index_best = cell_index

            prev_cell_indices += [cell_index_best]

        print(iter, x_best, calc_objective_func(x_best))

        x_preferred = x_best
        x_others = plane.get_vertices() + plane.get_center()

        optimizer.submit_data(x_preferred, x_others)

    x_opt = optimizer.retrieve_search_plane().get_center()
    f_opt = calc_objective_func(x_opt)

    print(x_opt, f_opt)
