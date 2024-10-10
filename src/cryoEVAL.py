import subprocess
import argparse
import re
import os


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
    # Get the full path of the script
    script_path = os.path.abspath(__file__) 
    # Extract the directory of the script
    script_dir = os.path.dirname(script_path)
    
    USALIGN = f'{script_dir}/utils/USalign'

    command = [
        USALIGN, 
        parsed_args.predicted_structure, 
        parsed_args.target_structure, 
        "-mm", "1", 
        "-ter", "0",
        "-d", "3.0",
        "-outfmt", "-1",  
        ]
    
    
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        output_log = result.stdout
                    
        if parsed_args.output_file is not None:
            with open(parsed_args.output_file, 'a') as f:
                f.write("*" * 52 + "\n")
                f.write("Credit from Zhang Lab: US-align (Version 20230609)\n")
                f.write("Reference: C Zhang, M Shine, AM Pyle, Y Zhang. (2022) Nat Methods\n")
                f.write("*" * 52 + "\n")
                f.write(output_log)


        ## Parse output ##
        aligned_length_pattern = r"Aligned length= (\d+)"
        tm_score_pattern = r"TM-score= ([0-9.]+) \(normalized by length of Structure_2"
        len_predict_pattern = r"Length of Structure_1: (\d+) residues"
        len_target_pattern = r"Length of Structure_2: (\d+) residues"
        
        # Search for 'Aligned length all d', not constrained by the distance threshold
        aligned_length_match = re.search(aligned_length_pattern, output_log)
        if aligned_length_match:
            aligned_length_alld = int(aligned_length_match.group(1))
        else:
            # Raise an exception if 'Aligned length' is not found
            raise ValueError("Aligned length not found in the output")

        # Search for 'TM-score'
        tm_score_match = re.search(tm_score_pattern, output_log)
        if tm_score_match:
            tm_score = float(tm_score_match.group(1))
        else:
            # Raise an exception if 'TM-score' is not found
            raise ValueError("TM-score not found in the output")
        
        # Search for 'Length of predicted and target structures'
        len_predict_match = re.search(len_predict_pattern, output_log)
        len_target_match = re.search(len_target_pattern, output_log)
        if len_predict_pattern and len_target_pattern:
            len_predict = int(len_predict_match.group(1))
            len_target = int(len_target_match.group(1))
        else:
            raise ValueError("Length of Predicted structure or Length of Target structure not found in the output")
        
        # extract the paired residues
        start_marker = '(":" denotes residue pairs of d < 3.0 Angstrom, "." denotes other aligned residues)'
        end_marker = "#Total CPU time is"
        start_index = output_log.find(start_marker) + len(start_marker)
        end_index = output_log.find(end_marker)
        extracted_text = output_log[start_index:end_index].strip()
        lines = extracted_text.split('\n')
        
        pred_seq = lines[0]
        target_seq = lines[2]
        match_notation = lines[1]
        
        assert len(pred_seq) == len(target_seq) == len(match_notation)
        assert match_notation.count(':') + match_notation.count('.') == aligned_length_alld
        
        aligned_length = match_notation.count(':')
        
        # calculate precision and recall and f1 score
        precision = aligned_length / len_predict
        recall = aligned_length / len_target
        f1score = 2 * precision * recall / (precision + recall)
        
        # count the number of matched residues that have the same type of amino acid
        aa_match = 0
        for i in range(len(pred_seq)):
            if (pred_seq[i] == target_seq[i]) and match_notation[i] == ':':
                aa_match += 1
                # print(pred_seq[i], target_seq[i], match_notation[i])
        
        # calculate sequence match and sequence recall       
        residue_match = aa_match / aligned_length
        residue_recall = aa_match / len_target
        
        # TMRR-score
        tmrr_score = 2 * tm_score * residue_recall / (tm_score + residue_recall)
        
        # completeness
        completeness = len_predict / len_target
        
        output = {
            "tm_score": "{:.3f}".format(tm_score),
            "aligned_length": aligned_length,
            "len_predict": len_predict,
            "len_target": len_target,
            "precision": "{:.3f}".format(precision),
            "recall": "{:.3f}".format(recall),
            "f1score": "{:.3f}".format(f1score),
            "residue_match": "{:.3f}".format(residue_match),
            "residue_recall": "{:.3f}".format(residue_recall),
            "tmrr_score": "{:.3f}".format(tmrr_score),
            "completeness": "{:.3f}".format(completeness),
        }
        
        if parsed_args.output_file is not None:
            with open(parsed_args.output_file, 'a') as f:
                f.write(f"\n")
                for key, value in output.items():
                    f.write(f"{key}: {value}\n")
        
        if parsed_args.verbose:
            print(output_log)
            for key, value in output.items():
                print(f"{key}: {value}")
                
        return output

    except subprocess.CalledProcessError as e:
        # Handle errors if the command fails
        print("Error:", e.stderr)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = add_args(parser)
    parsed_args = parser.parse_args()
    main(parsed_args)
