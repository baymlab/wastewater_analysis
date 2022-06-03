#!/usr/bin/env python3

import sys
import os
import copy
import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker
from mpl_toolkits import mplot3d
from matplotlib import cm
from math import log10, floor

superscript = str.maketrans("-0123456789.", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹·")

def main():
    parser = argparse.ArgumentParser(description="Evaluate predicted frequencies.")
    parser.add_argument('predictions', type=str, nargs='+', help="prediction files")
    parser.add_argument('--voc', dest='voc', type=str, required=True, help="comma-separated list of strains of interest")
    parser.add_argument('-o, --outdir', dest='outdir', required=True)
    parser.add_argument('--suffix', dest='suffix', default="", help="add suffix to output figure names")
    parser.add_argument('-v, --verbose', dest='verbose', action='store_true')
    parser.add_argument('--min_err', dest='min_err', default=0, type=float, help="minimal error (any samples with true error below this threshold are skipped; any predictions below this threshold are considered absent)")
    parser.add_argument('--min_ab', dest='min_ab', default=0, type=float, help="minimal abundance (any samples with true abundance below this threshold are skipped; any predictions below this threshold are considered absent)")
    parser.add_argument('--plot_error_value', dest='plot_err_val', default=1, type=float, help="error value for plots in which the error is not plotted. (Value is set to closest datapoint). Default: 1")
    parser.add_argument('--plot_abundance_value', dest='plot_ab_val', default=1, type=float, help="abundance value for plots in which the abundance is not plotted. (Value is set to closest datapoint). Default: 1")
    parser.add_argument('--no_plots', action='store_true')
    parser.add_argument('--output_format', dest='output_format', default='png', help="comma-separated list of desired output formats")
    parser.add_argument('--chimeric', action="store_true", help="indicate dataset with chimeric reads")
    parser.add_argument('--font_size', dest='font_size', default=12, type=int, help="set font size for the plots")
    args = parser.parse_args()

    false_pos_count = 0
    false_neg_count = 0
    true_pos_count = 0
    true_neg_count = 0
    err_list = []
    variant_set = set()
    voc_list = args.voc.split(',')
    output_formats = args.output_format.split(',')

    os.makedirs(args.outdir, exist_ok=True)

    # read predictions
    for filename in args.predictions:
        dir_name = filename.split('/')[-2]
        voc_name = dir_name.split('_')[0]
        err_freq = float(dir_name.split('_')[-1].lstrip('er' if not args.chimeric else 'ch'))
        voc_freq = float(dir_name.split('_')[-2].lstrip('ab'))
        if voc_name not in voc_list:
            continue
        elif err_freq < args.min_err:
            continue
        elif voc_freq < args.min_ab:
            continue
        variant_set.add(voc_name)
        with open(filename, 'r') as f:
            variant_found = False
            err_tups = []
            positives = []
            for line in f:
                if line[0] == "#":
                    continue
                [variant, tpm, ab, corrected_ab] = line.rstrip('\n').split('\t')
                if variant not in voc_list:
                    print("nah")
                    continue
                ab = float(ab)
                # Remove ↓ ?
                if ab < args.min_ab:
                    continue
                abs_err = abs(ab - voc_freq)
                positives.append(variant)
                if variant == voc_name:
                    variant_found = True
                    err_tups.append((voc_name, voc_freq, err_freq, abs_err, ab))
                else:
                    false_pos_count += 1
                    if args.verbose:
                        print("False positive: {} predicted at {}% with {}% error in {}".format(
                                variant, ab, err_freq, filename))
            if variant_found:
                true_pos_count += 1
                if len(err_tups) == 1:
                    err_list.append(err_tups[0])
                else:
                    voc_name = err_tups[0][0]
                    voc_freq = err_tups[0][1]
                    ab = sum([x[4] for x in err_tups])
                    abs_err = abs(ab - voc_freq)
                    err_list.append((voc_name, voc_freq, err_freq, abs_err, ab))
            else:
                false_neg_count += 1
                if args.verbose:
                    print("VOC not found in {}".format(filename))
                # add zero estimate to error list?
                # err_list.append((voc_name, voc_freq, err_freq, voc_freq, 0))
            for variant in voc_list:
                if variant not in positives and variant != voc_name:
                    true_neg_count += 1
            true_neg_count += len([x for x in voc_list if
                                    x not in positives and x != voc_name])


    # compute stats
    average_rel_err = sum([x[3]/x[1]*100 for x in err_list]) / len(err_list)
    average_rel_err_tp = (sum([x[3]/x[1]*100 for x in err_list if x[4] > 0])
                            / len(err_list))
    # print("average relative error: {}%".format(average_rel_err))
    print("average relative error of true positives: {}%".format(
                                                            average_rel_err_tp))
    print("total # true positives: {}".format(true_pos_count))
    print("total # true negatives: {}".format(true_neg_count))
    print("total # false positives: {}".format(false_pos_count))
    print("total # false negatives: {}".format(false_neg_count))

    fpr = save_dev(false_pos_count, (false_pos_count + true_neg_count))
    fnr = save_dev(false_neg_count, (false_neg_count + true_pos_count))
    recall = save_dev(true_pos_count, (true_pos_count + false_neg_count))
    precision = save_dev(true_pos_count, (true_pos_count + false_pos_count))
    print("FPR = {}".format(fpr))
    print("FNR = {}".format(fnr))
    print("Precision = {}".format(precision))
    print("Recall = {}\n".format(recall)) # sensitivity

    if args.no_plots:
        sys.exit()

    # sort error tuples by voc frequency
    err_list.sort(key = lambda x : x[1])
    err_list_err = copy.deepcopy(err_list)
    # sort second list by err frequency
    err_list_err.sort(key = lambda x : x[2])
    variant_list = sorted(list(variant_set))

    for voc in variant_list:
        os.makedirs(args.outdir + '/' + voc, exist_ok=True)

    # fix color per voc
    # print(len(variant_list))
    # colormap = cm.get_cmap('tab10', len(variant_list))
    # print(colormap)
    colors = {voc : cm.tab10((i))
                for i, voc in enumerate(variant_list)}

    # if args.joint_average:
    #     # compute average error for jointly evaluated VOCs
    #     if joint_voc_list != [""]:
    #         err_tups = [x for x in err_list if x[0] in joint_voc_list]
    #         new_err_list = [x for x in err_list if x[0] not in joint_voc_list]
    #         joint_voc_name = '/'.join(joint_voc_list)
    #         voc_freq = 0
    #         voc_freq_errs = []
    #         voc_freq_abs = []
    #         for tup in err_tups:
    #             if voc_freq > 0 and voc_freq != tup[1]:
    #                 av_err = sum(voc_freq_errs) / len(voc_freq_errs)
    #                 av_ab = sum(voc_freq_abs) / len(voc_freq_abs)
    #                 new_err_list.append((joint_voc_name, voc_freq, av_err, av_ab))
    #                 voc_freq_errs = []
    #                 voc_freq_abs = []
    #             voc_freq = tup[1]
    #             voc_freq_errs.append(tup[3])
    #             voc_freq_abs.append(tup[4])
    #         # finally add last tuple
    #         if voc_freq_errs:
    #             av_err = sum(voc_freq_errs) / len(voc_freq_errs)
    #             av_ab = sum(voc_freq_abs) / len(voc_freq_abs)
    #             new_err_list.append((joint_voc_name, voc_freq, av_err, av_ab))
    #         new_err_list.sort(key = lambda x : x[1])
    #         err_list = new_err_list
    #         variant_list = [voc for voc in variant_list if voc not in joint_voc_list]
    #         variant_list.append(joint_voc_name)
    #         colors[joint_voc_name] = colors[joint_voc_list[0]]

    _, f, e, _, _ = zip(*err_list_err)
    unique_err_vals = list(set(e))
    unique_freq_vals = list(set(f))

    # print(unique_freq_vals)
    # print("\n\n\n")
    # print(unique_err_vals)

    plot_single_err_val = min(unique_err_vals, key=lambda x:abs(x-args.plot_err_val))
    plot_single_frq_val = min(unique_freq_vals, key=lambda x:abs(x-args.plot_ab_val))


    plt.rcParams.update({'font.size': args.font_size}) # increase font size

    # abundance to relative error
    plt.figure()
    for voc in variant_list:
        freq_values = [x[1] for x in err_list if x[0] == voc and x[2] == plot_single_err_val]
        err_values = [x[3]/x[1]*100 for x in err_list if x[0] == voc and x[2] == plot_single_err_val]
        plt.plot(freq_values, err_values, label=voc, color=colors[voc])
        if (freq_values[0] > args.min_ab):
            plt.plot(freq_values[0], err_values[0], marker="s", color=colors[voc], markersize=6)
    # plt.vlines(plot_single_err_val, 0, 100)
    plt.legend()
    plt.grid(which="both", alpha=0.2)
    plt.ylim(-5, 105)
    plt.xlabel("True VoC frequency (%)")
    plt.ylabel("Relative prediction error (%)")
    # plt.gcf().set_size_inches(4, 3)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_error_plot{}.{}".format(args.outdir,
                                                     args.suffix,
                                                     format))

    # also plot on log scale
    plt.xscale('log')
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_error_plot_logscale{}.{}".format(args.outdir,
                                                              args.suffix,
                                                              format))


    # error to relative error
    plt.figure()
    for voc in variant_list:
        err_freq = [x[2] for x in err_list_err if x[0] == voc and x[1] == plot_single_frq_val]
        err_values = [x[3]/x[1]*100 for x in err_list_err if x[0] == voc and x[1] == plot_single_frq_val]
        plt.plot(err_freq, err_values, label=voc, color=colors[voc])
        if (err_freq[-1] < max(unique_err_vals)):
            plt.plot(err_freq[-1], err_values[-1], marker="s", color=colors[voc], markersize=6)
    plt.legend()
    plt.grid(which="both", alpha=0.2)
    plt.ylim(-5, 105)
    plt.xlabel("Error frequency (%)")
    plt.ylabel("Relative prediction error (%)")
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/error_error_plot{}.{}"
            .format(args.outdir, args.suffix, format))

    plt.xscale('log')
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/error_error_plot_logscale{}.{}"
            .format(args.outdir, args.suffix, format))

    # plot true vs estimated frequencies on a scatterplot
    plt.figure()
    for voc in variant_list:
        freq_values = [x[1] for x in err_list if x[0] == voc and x[2] == plot_single_err_val]
        est_values = [x[4] for x in err_list if x[0] == voc and x[2] == plot_single_err_val]
        plt.scatter(freq_values, est_values, label=voc, alpha=0.7,
                    color=colors[voc], s=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.07, 150)
    plt.ylim(0.07, 150)
    plt.plot([0, 100], [0, 100], 'k-', lw=0.75)
    plt.legend(prop={'size': 12}) #ncol=len(variants_list),
    plt.grid(which="both", alpha=0.2)
    plt.xlabel("True VOC frequency (%)")
    plt.ylabel("Estimated VOC frequency (%)")
    # # Hide the right and top spines
    # ax = plt.gca()
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_scatter_loglog{}.{}".format(args.outdir,
                                                         args.suffix,
                                                         format))
    

    # plot scatter with error gradient for every voc
    for voc in variant_list:
        plt.figure()

        freq_values = [x[1] for x in err_list_err if x[0] == voc]
        est_values = [x[4] for x in err_list_err if x[0] == voc]
        err_rate = [x[2] for x in err_list_err if x[0] == voc]

        plt.scatter(freq_values, est_values, alpha=0.7, 
                    c=err_rate, s=20, cmap='viridis', norm=matplotlib.colors.LogNorm())
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(0.07, 150)
        plt.ylim(0.07, 150)
        plt.plot([0, 100], [0, 100], 'k-', lw=0.75)
        # plt.legend(prop={'size': 12}) #ncol=len(variants_list),
        plt.colorbar(label='Induced error rate (%)')
        plt.grid(which="both", alpha=0.2)
        plt.xlabel("True {} frequency (%)".format(voc))
        plt.ylabel("Estimated {} frequency (%)".format(voc))
        # plt.title(voc)
        # # Hide the right and top spines
        # ax = plt.gca()
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        plt.tight_layout()
        for format in output_formats:
            plt.savefig("{}/{}/freq_error_scatter_loglog{}.{}".format(args.outdir,
                                                            voc,
                                                            args.suffix,
                                                            format))

    # plt.show()

    # plot 3d scatter graph
    for voc in variant_list:
        plt.figure()
        ax = plt.axes(projection='3d')

        freq_values = [x[1] for x in err_list if x[0] == voc]
        err_rate = [x[2] for x in err_list if x[0] == voc]
        err_values = [x[3]/x[1]*100 for x in err_list if x[0] == voc]

        # # cap relative error to 100
        # for i, val in enumerate(err_values):
        #     if val > 100:
        #         err_values[i] = 100

        X0 = np.array(err_rate)
        X0[X0 == 0] = 0.1
        X = np.log10(X0)

        Y0 = np.array(freq_values)
        Y0[Y0 == 0] = 0.1
        Y = np.log10(Y0)

        ax.scatter(X, Y, err_values, c=err_values, cmap='viridis', depthshade=False)
        ax.view_init(azim=135)

        ax.xaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))

        ax.set_xlabel("Induced error rate (%)")
        ax.set_ylabel("True {} frequency (%)".format(voc))
        ax.set_zlabel("Relative prediction error (%)")

        ax.set_zlim(0,100)

        # plt.title(voc)

        # plt.tight_layout()
        for format in output_formats:
            plt.savefig("{}/{}/freq_error_error_scatter_3d{}.{}".format(args.outdir, voc, args.suffix, format))

    if len(unique_err_vals) <= 1 or len(unique_freq_vals) <= 1:
        return

    # plot 3d graph
    for voc in variant_list:
        plt.figure()
        ax = plt.axes(projection='3d')

        freq_values = [x[1] for x in err_list if x[0] == voc]
        err_rate = [x[2] for x in err_list if x[0] == voc]
        err_values = [x[3]/x[1]*100 for x in err_list if x[0] == voc]

        # # cap relative error to 100
        # for i, val in enumerate(err_values):
        #     if val > 100:
        #         err_values[i] = 100

        X0 = np.array(err_rate)
        X0[X0 == 0] = 0.1
        X = np.log10(X0)

        Y0 = np.array(freq_values)
        Y0[Y0 == 0] = 0.1
        Y = np.log10(Y0)

        Z0 = np.array(err_values)
        Z0[Z0 == 0] = 0.1
        Z = np.log10(Z0)
        
        # ax.plot_trisurf(err_rate, freq_values, err_values)
        ax.plot_trisurf(X, Y, err_values, cmap='viridis', edgecolor='none')
        ax.view_init(azim=135)

        ax.xaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
        # ax.xaxis.set_major_locator(mticker.MaxNLocater(integer=True))

        ax.set_xlabel("Induced error rate (%)")
        ax.set_ylabel("True {} frequency (%)".format(voc))
        ax.set_zlabel("Relative prediction error (%)")

        ax.set_zlim(0,100)

        # plt.title(voc)

        # ax.set_xticks([-1,0,1,2,3])

        # plt.tight_layout()
        for format in output_formats:
            plt.savefig("{}/{}/freq_error_error_3d{}.{}".format(args.outdir, voc, args.suffix, format))

    # flat 3d plot
    # for voc in variant_list:
    #     plt.figure()

    #     freq_values = [x[1] for x in err_list if x[0] == voc]
    #     err_rate = [x[2] for x in err_list if x[0] == voc]
    #     err_values = [x[3]/x[1]*100 for x in err_list if x[0] == voc]

    #     # get space
    #     space_y = np.linspace(0, 9, num=10)
    #     space_x = np.logspace(-1, 2, num=100)
    #     # space_x = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]
    #     # space_x = np.concatenate((np.linspace(0, 1, num=33, endpoint=False), np.linspace(1, 10, num=33, endpoint=False), np.linspace(10, 100, num=33, endpoint=True)))
    #     # print(len(space_x), space_x)

    #     grid_y2, grid_x2 = np.mgrid[0:10,0:100]
    #     grid_x, grid_y = np.meshgrid(space_x, space_y)

    #     zz = griddata((freq_values, err_rate), err_values, (grid_x2, grid_y2), method='nearest')
    #     # print(zz)
    #     plt.plot(freq_values, err_rate, 'k.', ms=1)
    #     plt.pcolor(zz)

    #     plt.xscale('log')
    #     # plt.yscale('log')
    #     plt.xlim(0.1,100)
    #     plt.ylim(0,10)

    #     plt.xlabel("True VOC abundance (%)")
    #     plt.ylabel("Error frequency (%)")
    #     plt.colorbar(label="Relative prediction error (%)")
        
    #     for format in output_formats:
    #         plt.savefig("{}/{}/freq_error_error_flat{}.{}".format(args.outdir, voc, args.suffix, format))

    # for voc in variant_list:
    #     freq_values = [x[1] for x in err_list_err if x[0] == voc]
    #     err_rate = [x[2] for x in err_list_err if x[0] == voc]
    #     err_values = [x[3]/x[1]*100 for x in err_list_err if x[0] == voc]

    #     l1 = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0]
    #     l2 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0]
    #     print(len(l1), len(l2))
    #     l3 = [(x, y) for x in l2 for y in l1]
    #     print(len(l3))
    #     _, f, e, _, _ = zip(*err_list_err)

    #     l4 = [x for x in l3 if x not in zip(f, e)]
    #     print(l4, "|\n", list(zip(f, e)))

    #     # print(err_list_err)
    #     # print(len(list(zip(f, e))), list(zip(f, e)))

    #     # print(len(err_list_err), err_list_err)
    #     Z = np.array(err_values).reshape(len(set(err_rate)), len(set(freq_values)))
    #     x = freq_values
    #     y = err_rate

    #     fig, zx = plt.subplots()
    #     zx.pcolormesh(x, y, Z, shading='flat', vmin=min(Z), vmax=max(Z))
    #     X, Y = np.meshgrid(x, y)
    #     zx.plot(X.flat, Y.flat, 'o', color='m')

    #     for format in output_formats:
    #         plt.savefig("{}/{}/flat_test{}.{}".format(args.outdir, voc, args.suffix, format))
    # plt.show()
    return

def log_tick_formatter(val, pos=None):
    return "10" + str(round_sig(val, 1)).translate(superscript)

def round_sig(x, sig=2):
    if x == 0: return x
    return round(x, sig-int(floor(log10(abs(x))))-1)

def save_dev(a: int, b: int):
    if a == 0:
        return 0
    else:
        return a / b

if __name__ == "__main__":
    sys.exit(main())
