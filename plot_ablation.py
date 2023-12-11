import pandas as pd
import re,sys
import matplotlib.pyplot as plt

def remove_ci(accuracy):
    # Use a regular expression to match the first percentage value
    match = re.match(r'(\d+\.\d+)%', accuracy)

    if match:
        first_percentage = match.group(1)
        return float(first_percentage)
    else:
        print("No percentage found in the input string")


def get_accuracies(file_name):
    data_dict = {}
    # Open and read the text file
    with open(file_name, 'r') as file:
        for line in file:
            line_dict = {}
            results = line.strip().split(', ')
            
            for result in results:
                split = result.split(': ')
                param_name = split[0]
                param = split[1]
                line_dict[param_name] = param
            
            if "-train-" in line_dict["Setting"]:
                data_dict["Train"] = remove_ci(line_dict["Acc"])
            elif "-val-" in line_dict["Setting"]:
                data_dict["Val"] = remove_ci(line_dict["Acc"])
            elif "-test-" in line_dict["Setting"]:
                data_dict["Test"] = remove_ci(line_dict["Acc"])
            else:
                print("FIX")
    return data_dict

def main():
    if len(sys.argv) != 3:
        print("Usage: python plot_ablation.py <arg1> <arg2>")
        print("arg1: True -> n_shot | False -> n_way")
        print("arg2: method (eg. maml)")
        return

    # Extract and use the command-line arguments
    n_shot_bool = bool(sys.argv[1])
    method = sys.argv[2]

    all_results = {}
    if n_shot_bool:
        param = "n_shot"
        
    else:
        param = "n_way"

    samples = {}
    samples["n_shot"] = [1, 5, 10, 20, 50, 100] # Nshots
    samples["n_way"] = [1, 5, 10, 20, 35, 50]  # Nways

    for sample in samples[param]:
        file_name = "fewshotbench_v2/checkpoints/atacseq_adult_promoter_{}_{}{}/results.txt".format(method, sample, param[2:])
        result = get_accuracies(file_name)
        all_results[sample] = result

    df = pd.DataFrame(all_results).T
    df["Train"].plot()
    plt.xlabel(param)
    plt.ylabel('Accuracy')
    plt.title('Accuracy vs {}'.format(param))
    plt.savefig("plots/atacseq_adult_promoter_{}_{}{}.png".format(method, sample, param[2:])) 

if __name__ == "__main__":
    main()