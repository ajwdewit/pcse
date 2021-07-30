import argparse

from . import process_experiment_collections


def create_parser():
    parser = argparse.ArgumentParser(description='Run PCSE experiment collections', prog="python -m exp")
    parser.add_argument('--output_dir', dest='output_dir', action='store', default=None,
                        help='Directory to write the output figures from the experiments', required=False,
                        )
    return parser



def main():
    parser = create_parser()
    args = parser.parse_args()
    process_experiment_collections.main(args.output_dir)


if __name__ == "__main__":
    main()