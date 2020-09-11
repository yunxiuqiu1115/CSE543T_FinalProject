import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import us


def plot_win(gp_path, stan_path):
    stan_data = pd.read_csv(stan_path)
    gp_data = pd.read_csv(gp_path)

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


    idx = np.argsort(rep_win)
    rep_win= np.array(rep_win)[idx]
    dem_win = np.array(dem_win)[idx]
    states = np.array(states)[idx]
    print(states)

    ind = np.arange(len(dem_win))  
    width = 0.5
    p1 = plt.barh(ind, dem_win, width, color='blue', alpha=0.7)
    p2 = plt.barh(ind, rep_win, width, color='red', left=dem_win, alpha=0.7)

    to_abbr = us.states.mapping('name', 'abbr')
    labels = []
    for state in states:
        if state == "Georgia Special":
            labels.append('GA S')
        else:
            labels.append(to_abbr[state])

    plt.xlabel('Winning Probability')
    plt.title('Winning Probability by party at day 42')
    plt.yticks(ind, labels)
    plt.xticks(np.linspace(0, 1, num=11))
    # plt.legend((p1[0], p2[0]), ('Dem', 'Rep'))
    ax = plt.axes()
    ax.xaxis.set_ticks_position('none') 
    ax.yaxis.set_ticks_position('none') 
    plt.savefig("plots/wp42.pdf")
    plt.show()
    plt.close()


def plot_vote(gp_path, stan_path):
    stan_data = pd.read_csv(stan_path)
    gp_data = pd.read_csv(gp_path)
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

    idx = np.argsort(rep_vote)
    dem_vote = np.array(dem_vote)[idx]
    rep_vote = np.array(rep_vote)[idx]
    states = np.array(states)[idx]

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
    plt.title('Forecasted Vote by party at day 42')
    plt.yticks(ind, labels)
    plt.xticks(np.linspace(0, 1, num=11))
    # plt.legend((p1[0], p2[0]), ('Dem', 'Rep'))
    plt.savefig("plots/vote42.png")
    plt.show()
    plt.close()


def plot_vote_percentage(gp_path, stan_path):
    stan_data = pd.read_csv(stan_path)
    gp_data = pd.read_csv(gp_path)
    stan_data = stan_data.rename(columns={"median":"m"})

    dem_win = []
    rep_win = []
    dem_vote = []
    rep_vote = []
    dem_u = []
    rep_u = []
    dem_l = []
    rep_l = []
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
            vote = data[(data.state==state) & (data.candidate==candidate)].m.values[0]
            party = gp_data[(gp_data.cycle==2020) & (gp_data.state==state) & (gp_data.candidate==candidate)].party.values[0]
            u = data[(data.state==state) & (data.candidate==candidate)].upper95.values[0]
            l = data[(data.state==state) & (data.candidate==candidate)].lower95.values[0]
            if party==1: 
                rep_vote.append(vote)
                rep_metadata.append((state,candidate))
                rep_u.append(u)
                rep_l.append(l)
                rep_win.append(win)
            else:
                dem_vote.append(vote)
                dem_metadata.append((state,candidate))
                dem_u.append(u)
                dem_l.append(l)
                dem_win.append(win)

    idx = np.argsort(dem_win)
    dem_win = np.array(dem_win)[idx]
    rep_win = np.array(rep_win)[idx]
    dem_vote = np.array(dem_vote)[idx]
    rep_vote = np.array(rep_vote)[idx]
    dem_u = np.array(dem_u)[idx]
    rep_u= np.array(rep_u)[idx]
    dem_l = np.array(dem_l)[idx]
    rep_l = np.array(rep_l)[idx]
    states = np.array(states)[idx]

    ind = np.arange(len(dem_vote))  
    p1 = plt.errorbar(ind, dem_vote, dem_ci, marker='x', mfc='blue',
         mec='blue', ms=2, mew=2)
    p2 = plt.errorbar(ind, rep_vote, rep_ci, marker='o', mfc='red',
         mec='red', ms=2, mew=2)

    # p1 = plt.scatter(ind, dem_vote, color='blue', marker='x', alpha=0.8)
    # p2 = plt.scatter(ind, rep_vote, color='red', marker='o', alpha=0.8)

    to_abbr = us.states.mapping('name', 'abbr')
    labels = []
    for state in states:
        if state == "Georgia Special":
            labels.append('GA S')
        else:
            labels.append(to_abbr[state])

    plt.ylabel('Vote')
    plt.title('Forecasted Vote by party at day 42')
    plt.xticks(ind, labels)
    plt.yticks(np.linspace(0, 1, num=11))
    plt.legend((p1, p2), ('Dem', 'Rep'))
    plt.savefig("plots/vote42.png")
    plt.show()
    plt.close()


def main():
    gp_path = "results/LOOGP_2020day42_12.csv"
    stan_path = "results/stan_LOOGP_2020day42_12.csv"
    plot_win(gp_path, stan_path)
    # plot_win(gp_path, stan_path)

if __name__ == "__main__":
    main()
