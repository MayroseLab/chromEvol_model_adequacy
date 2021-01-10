from defs import *
from utils import *


def get_counts(filename):
    """
    reads the .counts_edit file and extracts the counts
    :param filename: counts file (original or simulated)
    :return: list of counts
    """
    with open(filename, "r") as counts_file:
        counts = []
        for line in counts_file:
            line = line.strip()
            if line.startswith('>'):  # taxa name
                continue
            else:
                if line=="x" or line=="X":  # discard counts with x
                    continue
                counts.append(int(line))
    return counts


def match_counts_to_tree(tree_file, out_dir):
    """
    Matches tree file to counts, in case of missing taxa or missing counts.
    :param tree_file: mlAncTree (containing counts in tips labels) in NEWICK format
    :param out_dir: where the processed tree and counts should be written to.
    :return:(1) tree_wo_counts without X taxa and without counts in their tip names
            (2) tree_with_counts without X taxa and with counts in the tips
    """
    t = Tree(tree_file, format=1)
    t = prune_x_from_tree(t, out_dir)
    produce_tree_without_counts(t, out_dir)
    remove_internal_labels(out_dir)


def prune_x_from_tree(t, out_dir):
    """
    Prune taxa with Xs as counts
    :param t: tree object
    :param out_dir: output directory
    :return: pruned tree (unpruned if not needed)
    """
    tips_to_prune = []
    all_tips = []
    for leaf in t:
        all_tips.append(leaf.name)
        name_with_x = re.search(".*\-X", leaf.name)
        if name_with_x:
            tips_to_prune.append(leaf.name)
    t.prune(list(set(all_tips) - set(tips_to_prune)))
    t.write(format=1, outfile=out_dir + tree_with_counts)
    return t


def produce_tree_without_counts(t, out_dir):
    """
    trim the counts digits from the tip labels
    :param t: tree object
    :param out_dir: output directory
    :return: tree without counts in tips labels
    """
    for leaf in t:
        name = re.search("(.*)\-[\d]", leaf.name)
        leaf.name = name.group(1)
    t.write(format=1, outfile=out_dir + tree_wo_counts)
    return t


def remove_internal_labels(out_dir):
    """
    remove internal nodes labels and add ":-1" at the end of the tree to make it rooted
    :param out_dir:
    :return: NA
    """
    with open(out_dir + tree_wo_counts, "r") as tree:
        t = tree.readline()
    # t = re.sub(";", ":-1;", t)  # remove ";" and add ":-1;"
    t = re.sub(r"_N\d+-\d+_", r"", t)  # replace all _N\d+-\d+_ with ""
    with open(out_dir + tree_wo_counts, "w+") as tree:
        tree.write(t)


def create_freq_file(tmp_line, freq_file):
    """
    create root frequency file
    :param tmp_line: the line in the results file that contains the root frequencies
    :param freq_file: root frequencies file to be written to
    :return:
    """
    text = tmp_line.split()
    with open(freq_file, "w+") as root_freq:
        root_freq.write("\n".join(text))


def get_rate_parameters(res, params, d1, freq_file):
    """
    get the rate parameters from the results file
    :param res: file handler of results file
    :param params: rate parameters names
    :param d1: dictionary to be updated
    :param freq_file: root frequencies file to be written to
    :return: updated parameters dictionary
    """
    for line in res:
        if line.startswith("#"):
            continue
        if line.startswith(params):
            tmp = re.search("(.*)\t(.*)", line)
            if tmp:
                key = tmp.group(1)
                val = int(tmp.group(2)) if key == "BASE_NUMBER" else float(tmp.group(2))
                d1[key] = val  # key = name of parameter, val = parameter's value
        else:
            create_freq_file(line, freq_file)
            break
    return d1


def get_params(res_file, freq_file):
    """
    parses results file, writes a root frequency file and creates parameters dictionary
    :param res_file: chromEvol results file (chromEvol.res(=)
    :param freq_file: root frequencies file to be written to
    :return: parameters dictionary
    """
    d1 = dict.fromkeys(["LOSS_CONST", "GAIN_CONST", "DUPL", "BASE_NUMBER_R", "BASE_NUMBER", "HALF_DUPL"], None)
    d2 = dict.fromkeys(["_lossConstR", "_gainConstR", "_duplConstR", "_baseNumberR", "_baseNumber", "_demiPloidyR"], None)
    mapping_dict = {"_lossConstR": "LOSS_CONST", "_gainConstR": "GAIN_CONST", "_duplConstR": "DUPL", "_baseNumberR": "BASE_NUMBER_R",  "_baseNumber": "BASE_NUMBER", "_demiPloidyR": "HALF_DUPL"}
    params = tuple(d1.keys())
    with open(res_file, "r") as res:
        d1 = get_rate_parameters(res, params, d1, freq_file)
    for key in d2:
        d2[key] = d1[mapping_dict[key]]
    is_demi = extract_line_from_file(res_file, "#Half_duplication rate is same as duplication rate")  # if CONST_RATE_DEMI parameter needs to be updated
    if is_demi:
        d2["_demiPloidyR"] = d2["_duplConstR"]
    d2["_simulationsTreeLength"] = extract_line_from_file(res_file, "#total tree", True, True)
    return d2


#####################################################################################################
#####################################################################################################
#                                                                                                   #
#                                 SIMULATIONS FUNCTIONS                                             #
#                                                                                                   #
#####################################################################################################
#####################################################################################################


def initialize_defaults(ma_output_dir, max_for_sim, tree_path, freq_file):
    """
    initialize several parameters in parameters dictionary, to be printed in the simulations parameters file
    :param ma_output_dir: where the simulations will be written to
    :param max_for_sim: current maximum allowed and initial maximum computed
    :param tree_path: phylogeny
    :param freq_file: frequency file path
    :return: parameters dictionary initialized with initial values
    """
    d = dict()
    d["_mainType"] = "mainSimulate"
    d["_outDir"] = ma_output_dir
    d["_treeFile"] = tree_path
    d["_freqFile"] = freq_file
    d["_simulationsJumpsStats"] = "expStats.txt"
    d["_simulationsIter"] = 1000
    if get_nsims() != 1000:
        d["_simulationsIter"] = get_nsims()
    d["_maxChrNumForSimulations"] = max_for_sim
    d["_branchMul"] = 1
    return d


def create_control_file(control_file, ma_output_dir, max_for_sim, tree_path, freq_file, params_dict, orig_counts):
    """
    create chromEvol parameters file for simulations, based on the results file
    :param control_file: name of file
    :param ma_output_dir: where the simulations will be written to
    :param max_for_sim: current maximum allowed and initial maximum computed
    :param tree_path: phylogeny
    :param freq_file: frequency file path
    :param params_dict: parameters dictionary
    :param orig_counts: original counts
    :return: NA
    """
    d = initialize_defaults(ma_output_dir, max_for_sim, tree_path, freq_file)
    if params_dict["_baseNumber"] is not None:
        # there will either be _maxBaseTransition after the second run. If not then the range is appropriate.
        params_dict["_maxBaseTransition"] = params_dict.get("_maxBaseTransition", max(max(orig_counts) - min(orig_counts), params_dict.get("_baseNumber")))
    d.update(params_dict)
    with open(control_file, "w+") as cf:
        for key in d.keys():
            if d[key] is not None:
                cf.write(key + " " + str(d[key]) + "\n")


def get_initial_max_allowed(orig_counts, results_file):
    """
    get the initial maximum chromosome number allowed from the results file and the current maximum allowed to be used in the first iteration
    :param orig_counts: original counts
    :param results_file: chromEvol results file
    :return: current maximum allowed and initial maximum computed
    """
    real_max = max(orig_counts)
    init_max_for_sim = max(real_max, min(real_max*10, 200))
    max_allowed = extract_line_from_file(results_file, "max chromosome allowed", True, True)
    return [max_allowed, init_max_for_sim]


def update_max_for_sim(m, init_max, max_allowed):
    """
    updates the current maximal number allowed by a factor
    :param m: multiplication factor
    :param init_max: initial maximum chromosome number allowed
    :param max_allowed: previous maximum chromosome number allowed
    :return: the current updated maximum for next iteration
    """
    max_for_sim = 200 * m + init_max
    if max_for_sim < max_allowed:
        max_for_sim = max_allowed
    return max_for_sim


def check_upper_bound(ma_output_dir, working_dir, max_for_sim, m):
    """
    checks if simulations reached upper bound of maximal chromosome number allowed
    :param ma_output_dir: where the simulations will be written to
    :param working_dir: output directory supplied by the user
    :param max_for_sim: current allowed maximal number
    :param m: multiplication factor
    :return: did simulations reach upper bound? True or False
    """
    sims_dirs = [x[0] for x in os.walk(ma_output_dir)]
    for i in range(1, len(sims_dirs)):
        sim_events_file = sims_dirs[i] + "/simEvents.txt"
        tmp = extract_line_from_file(sim_events_file, "Total number of transitions to max chromosome", True, True)
        if tmp:
            if tmp > 0:
                with open(working_dir + "/increasing_max_chr.txt", "w") as fh:
                    fh.write("Iteration number " + str(m+1) + ", max number is currently " + str(max_for_sim))
                return False
    return True


def run_simulations(results_file, ma_output_dir, orig_counts, tree_path, freq_file, params_dict, user_out_dir):
    """
    creates parameters file for simulations.
    checks upper bound of maximal chromosome number for simulations and updates accordingly.
    :param results_file: chromEvol results file
    :param ma_output_dir: where the simulations will be written to
    :param orig_counts: originla counts
    :param tree_path: phylogney
    :param freq_file: frequency file path
    :param params_dict: parameters dictionary
    :param user_out_dir: user output directory
    :return: NA
    """
    [max_allowed, init_max_for_sim] = get_initial_max_allowed(orig_counts, results_file)
    for mult in range(5):  # mult is the factor increasing _maxChrNumForSimulations
        max_for_sim = update_max_for_sim(mult, init_max_for_sim, max_allowed)
        create_control_file(ma_output_dir + sim_control, ma_output_dir, max_for_sim, tree_path, freq_file, params_dict, orig_counts)
        os.system('"' + chromevol_path + '" ' + ma_output_dir + sim_control)
        tmp = check_upper_bound(ma_output_dir, user_out_dir, max_for_sim, mult)
        if tmp:  # did not hit upper bound, no need to increase max_for_sim again
            break


#####################################################################################################
#####################################################################################################
#                                                                                                   #
#                                 BASE NUMBER SECOND RUN FUNCTIONS                                  #
#                                                                                                   #
#####################################################################################################
#####################################################################################################


def test_max_on_tree(base_num, counts_file, tree_file):
    """
    calculates the maximal transition on the inferred tree
    :param base_num: base number as parsed from the results file
    :param counts_file: counts input file
    :param tree_file: ml inferred tree
    :return: if the maximal transitions if equal or larger than the original range -
                return 0
            otherwise
                return a list of the new base number and the max base transition
    """
    counts = get_counts(counts_file)
    counts_range = range_of_lst(counts)
    max_base_on_tree = get_max_transition(tree_file)
    if max_base_on_tree >= counts_range:
        return 0
    max_base = max(max_base_on_tree, 3)
    base_num = min(max_base, base_num)
    return [base_num, max_base]


def create_control_file_second_run(control_file, out_dir, counts_file, tree_file, dupl_flag, bn, mb):
    """
    creates control file for second chromEvol run (for base_num models only)
    :param control_file: parameters file for chromEvol run
    :param out_dir: output directory of the second chromEvol run
    :param counts_file: counts file for the second chromEvol run
    :param tree_file: tree file for the second chromEvol run
    :param dupl_flag: extracted from the parameters dictionary - if None there is no duplication
    :param bn: base number
    :param mb: max base transition
    :return: NA
    """
    with open(control_file, "w+") as cf:
        cf.write("_mainType Optimize_Model\n")
        paths_params(cf, out_dir, counts_file, tree_file)
        fixed_params(cf)
        cf.write("_baseNumber " + str(bn) + "\n")
        cf.write("_maxBaseTransition " + str(mb) + "\n")
        if dupl_flag is not None:
            cf.write("_duplConstR 1\n")


def paths_params(control_file, out_dir, counts_file, tree_file):
    """
    Adds paths to the control file for the second chromEvol run
    :param control_file: file handler of control file
    :param out_dir: output directory for the second run
    :param counts_file: counts file for the second run
    :param tree_file: tree file for the second run
    :return: NA
    """
    control_file.write("_outDir " + out_dir + "\n")
    control_file.write("_dataFile " + counts_file + "\n")
    control_file.write("_treeFile " + tree_file + "\n")


def fixed_params(control_file):
    """
    Adds fixed parameters to the control file for the second chromEvol run
    :param control_file: file handler of control file
    :return: NA
    """
    control_file.write("_logFile log.txt\n")
    control_file.write("_maxChrNum -10\n")
    control_file.write("_minChrNum -1\n")
    control_file.write("_branchMul 999\n")
    control_file.write("_simulationsNum 1000\n")
    control_file.write("_logValue 6\n")
    control_file.write("_maxOptimizationIterations 5\n")
    control_file.write("_epsilonLLimprovement 0.01\n")
    control_file.write("_optimizePointsNum 10,2,1\n")
    control_file.write("_optimizeIterNum 0,1,3\n")
    control_file.write("_gainConstR 1\n")
    control_file.write("_lossConstR 1\n")
    control_file.write("_baseNumberR 1\n")
    control_file.write("_bOptBaseNumber 1\n")


def second_run(params_dict, results_path, counts_file, tree_file):
    """
    checks if a second chromEvol run is needed based on the maximal transition on the phylogeny
    :param params_dict: parameters dictionary
    :param results_path: where to print results to
    :param counts_file: original coutns file
    :param tree_file: phylogeny
    :return: NA
    """
    res = test_max_on_tree(params_dict.get("_baseNumber"), counts_file, results_path + tree_with_counts)
    if isinstance(res, list):  # need to re-run chromEvol with updated parameters
        create_control_file_second_run(results_path + "/second_run.params", results_path, counts_file, tree_file, params_dict.get("_duplConstR"), res[0], res[1])
        print("\nRunning another chromEvol optimization run\n")
        os.system('"' + chromevol_path + '" ' + results_path + "/second_run.params")
    open(results_path + "/second_run_tested", 'a').close()
