"""
Main function to run everything!
Discretize, extract, match, in that order.
"""
import lpath
from lpath.io import make_dir
from ._logger import Logger

log = Logger().get_logger(__name__)


def main(arguments):
    """
    Main function to run it all!

    Parameters
    ----------
    arguments : argparse.Namespace
        A namespace with all the parameter arguments

    """
    from lpath import discretize, extract, match, plot

    discretize.main(arguments)
    extract.main(arguments)
    match.main(arguments)
    plot.main(arguments)


def entry_point():
    """
    Entry point for discretize, extract, match steps

    """
    from lpath import argparser
    from lpath import discretize, extract, match, plot

    log.info(f'Running LPATH version {lpath.__version__}.')
    # Creating the subparsers for each subcommand
    subparsers = []
    parser = argparser.create_parser()
    parser, subparsers = lpath.argparser.create_subparsers(parser, subparsers)

    # Functions to be mapped to each subparser
    functions = [discretize.main, extract.main, match.main, plot.main, main]
    for subparser, func in zip(subparsers, functions):
        subparser.set_defaults(func=func)

    argparser.check_argv()

    # print(parser.__dict__)
    args = argparser.process_args(parser)
    log.info(f'LPATH arguments: {args}')

    make_dir(args)

    # print(args)
    # Run whatever function given
    args.func(args)
