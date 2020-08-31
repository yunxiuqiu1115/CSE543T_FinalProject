import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import us


def plot_win():
    stan_data = pd.read_csv("results/stan_sb20_all90_36.csv")
    gp_data = pd.read_csv("results/SB1992-2020all90_36.csv")

    dem_win = []
    rep_win = []
    dem_metadata = []
    rep_metadata = []
    states = []

    for state in stan_data.state.unique():
        data = stan_data[stan_data.state==state]
        if len(data.candidate.unique())>2:
            continue
        states.append(state)
        for candidate in data.candidate.unique():
            win = data[(data.state==state) & (data.candidate==candidate)].win.values[0]
            party = gp_data[(gp_data.cycle==2020) & (gp_data.state==state) & (gp_data.candidate==candidate)].party.values[0]
            if party==1: 
                rep_win.append(win)
                rep_metadata.append((state,candidate))
            else:
                dem_win.append(win)
                dem_metadata.append((state,candidate))

    ind = np.arange(len(dem_win))  
    width = 0.5
    p1 = plt.barh(ind, dem_win, width, color='blue', alpha=0.8)
    p2 = plt.barh(ind, rep_win, width, color='red', left=dem_win, alpha=0.8)

    to_abbr = us.states.mapping('name', 'abbr')
    labels = []
    for state in states:
        if state == "Georgia Special":
            labels.append('GA S')
        else:
            labels.append(to_abbr[state])

    plt.ylabel('Win')
    plt.title('Winning Probability by party at day 90')
    plt.yticks(ind, labels)
    plt.xticks(np.linspace(0, 1, num=11))
    plt.legend((p1[0], p2[0]), ('Dem', 'Rep'))
    plt.savefig("plots/wp90.png")
    plt.show()
    plt.close()

    for i in range(len(dem_win)):
        print(dem_metadata[i][0])
        print(dem_win[i])


def plot_vote():
    stan_data = pd.read_csv("results/stan_sb20_all90_36.csv")
    gp_data = pd.read_csv("results/SB1992-2020all90_36.csv")
    stan_data = stan_data.rename(columns={"median":"m"})

    dem_vote = []
    rep_vote = []
    dem_metadata = []
    rep_metadata = []
    states = []

    for state in stan_data.state.unique():
        data = stan_data[stan_data.state==state]
        if len(data.candidate.unique())>2:
            continue
        states.append(state)
        for candidate in data.candidate.unique():
            vote = data[(data.state==state) & (data.candidate==candidate)].m.values[0]
            party = gp_data[(gp_data.cycle==2020) & (gp_data.state==state) & (gp_data.candidate==candidate)].party.values[0]
            if party==1: 
                rep_vote.append(vote)
                rep_metadata.append((state,candidate))
            else:
                dem_vote.append(vote)
                dem_metadata.append((state,candidate))

    ind = np.arange(len(dem_vote))  
    width = 0.5
    p1 = plt.barh(ind, dem_vote, width, color='blue', alpha=0.8)
    p2 = plt.barh(ind, rep_vote, width, color='red', left=1-np.array(rep_vote), alpha=0.8)

    to_abbr = us.states.mapping('name', 'abbr')
    labels = []
    for state in states:
        if state == "Georgia Special":
            labels.append('GA S')
        else:
            labels.append(to_abbr[state])

    plt.ylabel('Vote')
    plt.title('Forecasted Vote by party at day 90')
    plt.yticks(ind, labels)
    plt.xticks(np.linspace(0, 1, num=11))
    plt.legend((p1[0], p2[0]), ('Dem', 'Rep'))
    plt.savefig("plots/vote90.png")
    plt.show()
    plt.close()

def main():
    plot_vote()
    plot_win()

if __name__ == "__main__":
    main()
