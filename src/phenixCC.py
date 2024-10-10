import subprocess
import argparse
import re

def add_args(parser):
    parser.add_argument(
        "--predicted-structure",
        "--p",
        "-p",
        type=str,
        required=True,
        help="Path to predicted structure",
    )
    parser.add_argument(
        "--target-structure",
        "--t",
        "-t",
        type=str,
        required=True,
        help="Path to target structure",
    )
    parser.add_argument(
        "--output-file",
        "--o",
        "-o",
        help="If set, path to save output file",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="If set, prints the results to the console",
    )
    
    return parser


def main(parsed_args):
    MatchScore = 'phenix.chain_comparison'
    
    command = [
        MatchScore, 
        parsed_args.target_structure, 
        parsed_args.predicted_structure, 
        # 'test_unique_part_of_target_only = False' 
        ]
    
    
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        output_log = result.stdout
        
        if parsed_args.output_file is not None:
            with open(parsed_args.output_file, 'a') as f:
                f.write("*" * 52 + "\n")
                f.write("Credit from PHENIX: phenix.chain_comparison function\n")
                f.write("*" * 52 + "\n")
                f.write(output_log)


        if parsed_args.verbose:
            print("*" * 52)
            print("Credit from PHENIX: phenix.chain_comparison function")
            print("*" * 52)
            print(output_log)
        
        ## Parse output ##
        pattern = r"_target\s+(\d+\.\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+)\s+(\d+)"

        match = re.search(pattern, output_log)
        
        if match:
            rmsd = float(match.group(1))
            close_n = int(match.group(2))
            far_n = int(match.group(3))
            forward = int(match.group(4))
            reverse = int(match.group(5))
            mixed = int(match.group(6))
            found = float(match.group(7))
            ca_score = float(match.group(8))
            seq_match = float(match.group(9))
            seq_score = float(match.group(10))
            mean_length = float(match.group(11))
            fragments = int(match.group(12))
            bad_connections = int(match.group(13))

            # print(f"RMSD: {rmsd}, CLOSE N: {close_n}, FAR N: {far_n}, FORWARD: {forward}, REVERSE: {reverse}, MIXED: {mixed}, CA Match: {ca_match}, CA SCORE: {ca_score}, SEQ MATCH: {seq_match}, SEQ SCORE: {seq_score}, MEAN LENGTH: {mean_length}, FRAGMENTS: {fragments}, BAD CONNECTIONS: {bad_connections}")
        else:
            raise ValueError("MatchScore output not found")
        
        
        output = {
            "found": found,
            "ca_score": ca_score,
            "seq_match": seq_match,
            "seq_score": seq_score,
            "mean_length":mean_length,
            "fragments": fragments,
            "bad_connections": bad_connections, 
        }

        return output
    
    except subprocess.CalledProcessError as e:
        # Handle errors if the command fails
        print("Error:", e.stderr)
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = add_args(parser)
    parsed_args = parser.parse_args()
    main(parsed_args)
