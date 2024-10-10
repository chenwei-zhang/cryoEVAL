import argparse
import os

from modelangeloEval import main as modelangeloEval_main
from phenixCC import main as phenixCC_main
from cryoEVAL import main as cryoEVAL_main


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
        required=True,
        help="If set, path to save output file",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="If set, prints the results to the console",
    )
    parser.add_argument(
        "--cryoEVAL",
        type=bool,
        default=True,
        help="If True, do cryoEVAL evaluation",
    )    
    parser.add_argument(
        "--modelangelo",
        type=bool,
        default=False,
        help="If True, do modelangelo evaluation",
    )
    parser.add_argument(
        "--phenix",
        type=bool,
        default=False,
        help="If True, do phenix.chain_comparison evaluation",
    )
    parser.add_argument(
        "--max-dist",
        type=float,
        default=3,
        help="In angstrom (Å), the maximum distance for correspondence",
    )
    parser.add_argument(
        "--match-type",
        default="both",
        choices=["both", "protein", "nucleotide"],
    )
    parser.add_argument(
        "--output-structure",
        help="If set, saves the sequence recall results to an mmCIF file, "
        "B-factors of 100 correspond to correct classifications and "
        "B-factors of 0 correspond to wrong classifications",
    )
    
    return parser


def main(parsed_args):
    
    given_output_file = parsed_args.output_file
    dir_path = os.path.dirname(given_output_file)
    file_name_without_ext = os.path.basename(given_output_file).split('.')[0]
        
    #### Write log file ####
    log_file = os.path.join(dir_path, f"{file_name_without_ext}_TraceLog.log")
    parsed_args.output_file = log_file    
    
    cryoEVAL_output = {}
    modelangelo_output = {}
    phenix_output = {}

    if parsed_args.cryoEVAL:
        print("Run cryoEVAL ...")
        with open(log_file, 'a') as f:
            f.write("[Measure - cryoEVAL] \n")
            f.flush()
            cryoEVAL_output = cryoEVAL_main(parsed_args)
            f.write("\nDONE!\n\n\n\n\n")
            f.flush()
    
    if parsed_args.modelangelo:
        print("Run ModelAngelo ...")
        with open(log_file, 'a') as f:
            f.write("[Measure - ModelAngelo]\n")
            f.flush()
            modelangelo_output = modelangeloEval_main(parsed_args)
            f.write("\nDONE!\n\n\n\n\n")
            f.flush()
            
    if parsed_args.phenix:
        print("Run PHENIX ...")
        with open(log_file, 'a') as f:
            f.write("[Measure - PHENIX] \n")
            f.flush()
            phenix_output = phenixCC_main(parsed_args)
            f.write("\nDONE!\n\n\n\n\n")
            f.flush()
    
    
    #### Write final output file ####
    parsed_args.output_file = given_output_file
        
    print("\nEvaluation Results:")
    
    with open(parsed_args.output_file, 'a') as f:
        
        # Interate over the dictionary and write key-value ot the file
        pred_name = os.path.basename(parsed_args.predicted_structure)
        f.write(f"# Prediction file: {pred_name} #\n")
        
        if parsed_args.cryoEVAL:
            # From cryoEVAL
            f.write("########## cryoEVAL Results #############\n")
            f.flush()
            for key, value in cryoEVAL_output.items():
                f.write(f"{key}: {value}\n") 
                print(f"  - {key}: {value}")  
        
        if parsed_args.modelangelo:
            # From ModelAngelo
            f.write("########## ModelAngelo Results #############\n")
            f.flush()
            for key, value in modelangelo_output.items():
                f.write(f"{key}: {value}\n") 
                print(f"  - {key}: {value}")
        
        if parsed_args.phenix:
            # From PHENIX
            f.write("########## PHENIX Results #############\n")
            f.flush()
            for key, value in phenix_output.items():
                f.write(f"{key}: {value}\n") 
                print(f"  - {key}: {value}")
                
        f.write(f"\n") 
        f.write(f"\n")
        f.flush()
        
    return (cryoEVAL_output, modelangelo_output, phenix_output)
            
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = add_args(parser)
    parsed_args = parser.parse_args()
    main(parsed_args)
    
    
    