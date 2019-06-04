import argparse
import sys


def print_includes():
    from libsequence import get_includes

    idirs = [get_includes()]

    print(' '.join('-I' + idir for idir in idirs))


def get_includes():
    from libsequence import get_includes
    print(get_includes())


def main():
    parser = argparse.ArgumentParser(prog='python -m libsequence')
    parser.add_argument('--includes', action='store_true',
                        help='Print the CPPFLAGS required to use libsequence headers.')
    parser.add_argument('--get_includes', action='store_true',
                        help='Print the location of the libsequence headers.')
    args = parser.parse_args()
    if not sys.argv[1:]:
        parser.print_help()
    if args.includes:
        print_includes()
        return
    if args.get_includes:
        get_includes()
        return


if __name__ == "__main__":
    main()
