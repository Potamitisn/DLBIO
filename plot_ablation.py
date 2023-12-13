import re, os
import matplotlib.pyplot as plt

def get_num_param(input_string):
    """
    Given a string of the form "Xword", where X is a number and word is a string,
    it splits the string into the number and the string.
    In our case the word/string is the parameter we are varying.
    """
    match = re.match(r'(\d+)([a-zA-Z]+)', input_string)

    if match:
        num = int(match.group(1))
        param = "n_" + match.group(2)
        return num, param
    else:
        print("No match found.")

def get_accuracies(file_name):
    """
    Reads the train, val and test accuracies from the results.txt file,
    as well as their confidence intervals.
    """
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
                data_dict["Train"] = float(line_dict["Acc"].split(" +- ")[0].rstrip("%"))
                data_dict["TrainCI"] = float(line_dict["Acc"].split(" +- ")[1].rstrip("%"))
            elif "-val-" in line_dict["Setting"]:
                data_dict["Val"] = float(line_dict["Acc"].split(" +- ")[0].rstrip("%"))
                data_dict["ValCI"] = float(line_dict["Acc"].split(" +- ")[1].rstrip("%"))
            elif "-test-" in line_dict["Setting"]:
                data_dict["Test"] = float(line_dict["Acc"].split(" +- ")[0].rstrip("%"))
                data_dict["TestCI"] = float(line_dict["Acc"].split(" +- ")[1].rstrip("%"))
            else:
                print("FIX")
    return data_dict

    

def main():
    """
    Plots the results of the ablation study.
    """
    results = {"n_shot": {}, "n_way": {}}
    checkpoints_dir = "fewshotbench_v2/checkpoints"
    # Get all the experiments that we run
    experiments = [f for f in os.listdir(checkpoints_dir) if os.path.isdir(os.path.join(checkpoints_dir, f))]

    # Get the results for each experiment
    for experiment in experiments:
        # Change name for readability
        experiment_name = experiment.replace("baseline_pp", "baseline++")
        method = experiment_name.split('_')[3]
        num, param = get_num_param(experiment_name.split('_')[-1])
        file_name = os.path.join(checkpoints_dir, experiment, "results.txt")
        accuracies = get_accuracies(file_name)
        if method not in results[param]:
            results[param][method] = [(num, accuracies["Test"], accuracies["TestCI"])]
        else:
            results[param][method].append((num, accuracies["Test"], accuracies["TestCI"]))

    # Plot the results
    fig, axes = plt.subplots(1, 2, figsize=(8,4), sharey=True)

    # Sort the results by the number of shots/ways and plot Test accuracy + CI
    for method, exp_results in results["n_way"].items():
        exp_results.sort(key=lambda x: x[0])
        x, y, ci = zip(*exp_results)
        lower_ci = tuple(y[i]-ci[i] for i in range(len(y)))
        upper_ci = tuple(y[i]+ci[i] for i in range(len(y)))
        axes[0].plot(x, y, label=method, linewidth=2, marker='o', markersize=5)
        axes[0].fill_between(x, lower_ci, upper_ci, alpha=0.3)
    
    for method, exp_results in results["n_shot"].items():
        exp_results.sort(key=lambda x: x[0])
        x, y, ci = zip(*exp_results)
        lower_ci = tuple(y[i]-ci[i] for i in range(len(y)))
        upper_ci = tuple(y[i]+ci[i] for i in range(len(y)))
        axes[1].plot(x, y, label=method, linewidth=2, marker='o', markersize=5)
        axes[1].fill_between(x, lower_ci, upper_ci, alpha=0.3)

    axes[1].legend()
    axes[0].legend()
    axes[1].set_xlabel("Shot")
    axes[0].set_xlabel("Way")
    axes[0].set_ylabel("Test accuracy")
    
    fig.suptitle('Ablation study', fontsize=14)
    plt.savefig("plots/waysVshots.png") 
    

if __name__ == "__main__":
    main()