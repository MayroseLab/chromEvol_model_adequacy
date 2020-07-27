from utils import *
from defs import *
from data_processing import *
from analysis import *

if __name__ == '__main__':
    counts_file, tree_file, results_file, user_output_dir = get_arguments()
    model_adequacy_out_dir = user_output_dir + "/adequacy_test/"
    if not os.path.exists(model_adequacy_out_dir):
        os.system("mkdir -p " + model_adequacy_out_dir)
    parameters = get_params(results_file, user_output_dir + root_freq_filename)
    match_counts_to_tree(tree_file, user_output_dir)
    original_counts = get_counts(counts_file)
    original_counts_statistics = calculate_statistics(original_counts, user_output_dir + tree_with_counts, model_adequacy_out_dir + "orig_stats", True)
    if parameters["_baseNumber"] is not None:  # base_num flag - execute second chromEvol run, if needed
        if not os.path.exists(user_output_dir + "/second_run_tested"):
            second_run(parameters, user_output_dir, counts_file, user_output_dir + tree_wo_counts)
    run_simulations(results_file, model_adequacy_out_dir, original_counts, user_output_dir + tree_wo_counts, user_output_dir + root_freq_filename, parameters, user_output_dir)
    model_adequacy(model_adequacy_out_dir, original_counts_statistics, user_output_dir, original_counts_statistics)
