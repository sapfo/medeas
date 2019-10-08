import numpy as np
from scipy.optimize import least_squares


def compute_coalescence_time(ts: 'np.array',ns: 'np.array')->'float':
    """Given a set of split time ts and the effective sizes ns, compute the expected
    Coalescence time. Notice that the length of ns should be 1 + lenght of ts"""
    ts = np.concatenate(([0], ts, [np.inf]))

    delta_ts = (ts[1:] - ts[:-1])
    nb_ns = len(ns)

    probability_unnormalized = 1. - np.exp(-delta_ts/ns)
    probability = np.full(nb_ns,0.)
    for index_n in range(nb_ns):
        probability_previous = np.sum(probability[0:index_n])
        new_prob = (1-probability_previous)*probability_unnormalized[index_n]
        probability[index_n] = new_prob

    all_contribution = (ts[:-1] + ns + delta_ts/(1-np.exp(delta_ts/ns)))*probability
    all_contribution[-1] = (ts[-2] + ns[-1])*probability[-1]
    coalescence_time = np.sum(all_contribution)
    return coalescence_time


def extrapolate_split_time(tree,
                           index_split_matrix,
                           all_coalescence_time_within,
                           all_coalescece_time_between
                           ):
    all_split_time = np.full(len(all_coalescece_time_between),0.)
    all_effective_size = np.full(len(all_coalescence_time_within), np.max(all_coalescence_time_within))

    Ds_above = np.array([])
    Ns_above = np.array([np.max(all_coalescence_time_within)])

    def set_sub_tree(tree, Ds_above, Ns_above):
        """Compute the effective size and the split time of a given split and all his sup split given the past history
         defined by the previous change in population size"""
        if tree.is_tip():
            return None
        left_tree = tree.children[0]
        right_tree = tree.children[1]

        def get_largest_id_and_value_tree(tree):
            largest_value = 0
            largest_id = 0
            if tree.is_tip():
                largest_id = int(tree.name)
                largest_value = all_coalescence_time_within[largest_id]
            else:
                for tips in tree.tips():
                    id = int(tips.name)
                    time = all_coalescence_time_within[id]
                    if time > largest_value:
                        largest_id = id
                        largest_value = time
            return (largest_id, largest_value)

        largest_left_id, largest_left_value = get_largest_id_and_value_tree(left_tree)
        largest_right_id, largest_right_value = get_largest_id_and_value_tree(right_tree)
        all_tree = [left_tree, right_tree]

        def compute_t_between_from_ds_and_ns(Dlargesmall, Nsmall) -> 'list':
            tlargelarge = compute_coalescence_time(Ds_above, Ns_above)
            tsmallsmall = compute_coalescence_time(np.append([Dlargesmall], Ds_above),
                                                   np.append([Nsmall], Ns_above))
            tlargesmall = Dlargesmall + compute_coalescence_time(Ds_above - Dlargesmall, Ns_above)
            return (tlargelarge, tsmallsmall, tlargesmall)

        all_largest_id = [largest_left_id, largest_right_id]
        all_largest_value = [largest_left_value, largest_right_value]
        which_largest_value = np.argmax([all_largest_value])
        smallest_id = all_largest_id[1-which_largest_value]
        tlargelarge = all_largest_value[which_largest_value]
        tsmallsmall = all_largest_value[1 - which_largest_value]
        tlargesmall = all_coalescece_time_between[index_split_matrix[largest_left_id, largest_right_id]]

        real_ts = (tlargelarge, tsmallsmall, tlargesmall)

        def dev(x) -> 'np.ndarray[float]':
            """Find deviation (residuals)"""
            putative_ts = np.array(compute_t_between_from_ds_and_ns(x[0], x[1]))
            return np.real(putative_ts - real_ts)

        res = least_squares(dev, (tlargesmall - tlargelarge, Ns_above[0]), gtol=1.e-15)

        all_split_time[index_split_matrix[largest_left_id, largest_right_id]] = res.x[0]
        all_effective_size[smallest_id] = res.x[1]

        set_sub_tree(all_tree[which_largest_value],
                     Ds_above,
                     Ns_above)
        set_sub_tree(all_tree[1 - which_largest_value],
                     np.insert(Ds_above, 0, res.x[0]),
                     np.insert(Ns_above, 0, res.x[1]))
    set_sub_tree(tree,
                 Ds_above,
                 Ns_above)
    return((all_effective_size, all_split_time))