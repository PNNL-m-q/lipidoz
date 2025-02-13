"""
lipidoz/_cli.py
Dylan Ross (dylan.ross@pnnl.gov)

    Command-line interface for running analysis
"""


import argparse

from lipidoz.workflows import (
    run_isotope_scoring_workflow, 
    save_isotope_scoring_results,
    write_isotope_scoring_report_xlsx
)
from lipidoz.gui.app import LozApp


def _setup_process_subparser(subparser: argparse.ArgumentParser):
    """ set up the subparser for process subcommand """
    subparser.add_argument(
        "OZ_DATA",
        help="OzID data file (.uimf or .mza)"
    )
    subparser.add_argument(
        "TARGET_LIST",
        help="target list file  (.csv)"
    )
    subparser.add_argument(
        "RESULTS_LOZ",
        help="save LipdOz results to file  (.loz)"
    )
    subparser.add_argument(
        "--results-xlsx",
        type=str,
        default=None,
        help="export LipdOz results to spreadsheet  (.xlsx)"
    )
    subparser.add_argument(
        "--mz-tol",
        type=float,
        default=0.05,
        help="m/z tolerance for chromatogram extraction (default=0.05)"
    )
    subparser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="print extra debugging messages"
    )
    subparser.add_argument(
        "--D-label",
        dest="d_lab",
        type=int,
        default=None,
        help="deuterium label"
    )
    subparser.add_argument(
        "--D-label-in-NL",
        dest="d_lab_nl",
        action="store_true",
        default=False,
        help="include deuterium label in neutral loss"
    )


def _setup_gui_subparser(subparser: argparse.ArgumentParser):
    """ set up the subparser for gui subcommand """
    subparser.add_argument(
        "--results-loz",
        type=str,
        default=None,
        help="view LipidOz results file (.loz)"
    )


def _setup_arg_parser():
    """ set up the argument parser """
    # set up main parser
    parser = argparse.ArgumentParser(prog="LipidOz")
    _subparsers = parser.add_subparsers(
        title="subcommands", 
        required=True,
        dest="subcommand"
    )
    # set up processing subparser
    _setup_process_subparser(
            _subparsers.add_parser(
            "process", 
            help="process OzID data"
        )
    )
    # set up gui subparser
    _setup_gui_subparser(
        _subparsers.add_parser(
            "gui", 
            help="start graphical user interface"
        )
    )
    return parser


def run():
    args = _setup_arg_parser().parse_args()
    match args.subcommand:
        case "process":
            results = run_isotope_scoring_workflow(
                args.OZ_DATA, 
                args.TARGET_LIST, 
                args.mz_tol,
                debug_flag="text" if args.debug else None,
                info_cb=print
            )
            save_isotope_scoring_results(results, args.RESULTS_LOZ)
            if args.results_xlsx is not None:
                write_isotope_scoring_report_xlsx(results, args.results_xlsx)
        case "gui":
            app = LozApp()
            if args.results_loz is not None:
                app._load_and_display_results(args.results_loz)
            else:
                app.run()