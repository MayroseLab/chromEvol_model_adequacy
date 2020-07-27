import os
import regex as re
import argparse
from ete3 import Tree

#  fixed names
log = "/log.txt"
CE_res_filename = "/chromEvol.res"
mlAncTree = "/mlAncestors.tree"
posterior_tree = "/posteriorAncestors.tree"
exp_tree = "/exp.tree"
expectation_file = "/expectations.txt"
anc_prob = "/ancestorsProbs.txt"
tree_with_counts = "/tree_with_counts.tree"
tree_wo_counts = "/tree_wo_counts.tree"
root_freq_filename = "/root_freq"
sim_control = "/param_sim"
adequacy_vec = "/adequacy_vec"


chromevol_path, nsims = None, None


def get_chromevol_path():
	return chromevol_path


def set_chromevol_path(path):
	global chromevol_path
	chromevol_path = path


def get_nsims():
	return nsims


def set_nsims(n):
	global nsims
	nsims = n


def get_arguments():
	parser = argparse.ArgumentParser(description='Pipeline to test model adequacy of a chromEvol model')
	parser.add_argument('--counts', '-c', help='Counts file', required=True)
	parser.add_argument('--tree', '-t', help='Tree file', required=True)
	parser.add_argument('--results_file', '-r', help='chromEvol results file', required=True)
	parser.add_argument('--user_output_dir', '-out', help='Output directory', required=True)
	parser.add_argument('--sims_per_tree', '-n', help='Number of simulations for the adequacy test', required=False, default=1000)
	parser.add_argument('--chromevol', '-ce', help='chromEvol executable path', required=True)

	# parse arguments
	args = parser.parse_args()
	counts_file = args.counts
	tree_file = args.tree
	results_file = args.results_file
	user_output_dir = args.user_output_dir
	set_nsims(int(args.sims_per_tree))
	set_chromevol_path(args.chromevol)

	return counts_file, tree_file, results_file, user_output_dir
