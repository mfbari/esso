import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file')
    parser.add_argument('output_file')
    args = parser.parse_args()

    with open(args.input_file) as in_file:
        file_content = in_file.read()
        with open(args.output_file, 'w') as out_file:
            out_file.write(file_content.replace(',', ' '))


if __name__ == '__main__':
    main()
