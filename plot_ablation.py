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


def get_df(param, method):
    all_results = {}


    samples = {}
    samples["n_shot"] = [1, 10, 20, 50, 100] # Nshots
    samples["n_way"] = [1, 10, 20, 35, 50]  # Nways

    for sample in samples[param]:
        file_name = "fewshotbench_v2/checkpoints/atacseq_adult_promoter_{}_{}{}/results.txt".format(method, sample, param[2:])
        result = get_accuracies(file_name)
        all_results[sample] = result

    df = pd.DataFrame(all_results).T
    return df


def main():
    methods = ["matchingnet", "protonet"]

    fig, axes = plt.subplots(1, 2, figsize=(8,4), sharey=True)
    
    for method in methods+["maml"]:
        df = get_df("n_shot", method)
        df["Test"].plot(ax=axes[1], label=method)

    for method in methods:
        df = get_df("n_way", method)
        df["Test"].plot(ax=axes[0], label=method)
    
    axes[1].legend()
    axes[0].legend()
    axes[1].set_xlabel("Shot")
    axes[0].set_xlabel("Way")
    axes[0].set_ylabel("Test accuracy")
    
    fig.suptitle('Ablation study', fontsize=14)
    plt.savefig("plots/waysVshots.png") 
    

if __name__ == "__main__":
    main()