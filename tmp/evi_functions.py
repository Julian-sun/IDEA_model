import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def func(m, series, r_a, mu, c):
    """This function calculates the EVIt-1,t which can be further used for the
    creation of the EVI ind. The imput funcion is: m: the rolling window size;
    serie: the time series where we desire to employ our result;
    r_a: the lag for the mooving average to be applied to the series;
    mu: the window for posterior mean (used in the original method as 7)
    c1 and c2 are treshold limits
    """

    tseries = pd.Series(series).rolling(window=r_a).mean()
    std_windows = []
    for i in range(0, len(tseries) - m + 1):
        win = tseries[0+i:m+i]
        std_windows.append(win.std())

    evi_t1_t = round(
        (pd.Series(std_windows) - pd.Series(std_windows).shift()) /
        abs(pd.Series(std_windows).shift()), 2
    )

    evi_t1_t_2 = [0] * (m-1) + evi_t1_t.to_list()

    mean_windows = []

    for i in range(0, len(tseries) - mu + 1):
        win = series[0+i:mu+i]

        mean_windows.append(win.mean())

    mu_ = [0] * (mu-1) + mean_windows

    # Create a dataframe with the results
    d = {
        'series': series,
        'tseries': tseries,
        'evi_t1_t': evi_t1_t_2,
        'mu': mu_
    }

    data = pd.DataFrame(data=d)

    data.evi_t1_t = data.evi_t1_t.fillna(0)

    data.mu = round(data.mu, 0)

    data.loc[(data['series'] < 0.75 * data['mu']), 'cat_mu'] = 0
    data.loc[((data['series'] >= 0.75 * data['mu']) &
              (data['series'] < data['mu'])), 'cat_mu'] = 0.5

    data.loc[(data['series'] >= data['mu']), 'cat_mu'] = 1

    data.loc[(data['evi_t1_t'] < c), 'cat_evi'] = 0
    # data.loc[(data['evi_t1_t'] >= c1) & (data['evi_t1_t'] < c2),'cat_evi'] = 0.5
    data.loc[(data['evi_t1_t'] >= c), 'cat_evi'] = 1

    data.loc[(data['cat_evi'] == 0), 'ind'] = 0
    data.loc[(data['cat_evi'] == 1) & (data['cat_mu'] == 0), 'ind'] = 0
    data.loc[(data['cat_evi'] == 1) & (data['cat_mu'] == 0.5), 'ind'] = 0
    # data.loc[(data['cat_evi'] == 0) & (data['cat_mu'] == 1),'ind'] = 0.5
    data.loc[(data['cat_evi'] == 1) & (data['cat_mu'] == 1), 'ind'] = 1

    data.loc[(data['cat_mu'] == 0) | (data['evi_t1_t'] < 0), 'ind_2'] = 0
    # data.loc[(data['evi_t1_t'] = 0) & (data['cat_mu'] == 1),'ind_2'] = 0
    data.loc[(data['evi_t1_t'] >= 0) &
             (data['evi_t1_t'] <= c) & (data['cat_mu'] >= 0.5), 'ind_2'] = 0.5

    data.loc[(data['evi_t1_t'] > c) & (data['cat_mu'] == 0.5), 'ind_2'] = 0.5
    data.loc[(data['evi_t1_t'] > c) & (data['cat_mu'] == 1), 'ind_2'] = 1

    return data


def plot_func(dtf, c):

    x = list(range(0, len(dtf)))

    dtf = dtf.assign(Days=x)

    fig = plt.figure(figsize=(15, 3))
    plt.plot(x, dtf.series)
    plt.plot(x, dtf.tseries, '-o')
    plt.plot(x, dtf.mu, linestyle='--', label='mu')
    # plt.plot(x, 0.75*dtf.mu,linestyle = '--',label = '0.75mu')
    plt.legend()
    plt.xticks(rotation=90)
    plt.tight_layout()

    # plot2

    fig = plt.figure(figsize=(15, 3))
    plt.plot(x, dtf.evi_t1_t)
    plt.axhline(y=c, color='y', linestyle='--', label='c')
    # plt.axhline(y = c2, color = 'r', linestyle = '--', label='c2')
    plt.legend()
    plt.ylim([-1, 1])
    plt.xlabel('Days')
    plt.ylabel('EVI_t-1,t')
    plt.xticks(rotation=90)
    plt.tight_layout()

    # plot3
    fig = plt.figure(figsize=(15, 4))
    sns.scatterplot(x="Days", y="tseries", data=dtf, hue="ind",
                    palette="RdYlBu_r")  # gist_gray, flare,RdYlBu_r

    # plot4
    fig = plt.figure(figsize=(15, 4))
    sns.scatterplot(x="Days", y="tseries", data=dtf, hue="ind_2",
                    palette="RdYlBu_r")  # gist_gray, flare,RdYlBu_r
    plt.plot(x, dtf.mu, linestyle='--', label='mu')


def plot_func2(dtf, c):

    # x = list(range(0, len(dtf)))
    x = dtf.year_week_ts
    # dtf = dtf.assign(Days = x)

    fig = plt.figure(figsize=(15, 3))
    plt.plot(x, dtf.atend_ivas)
    plt.plot(x, dtf.tseries, '-o')
    plt.plot(x, dtf.mu, linestyle='--', label='mu')
    # plt.plot(x, 0.75*dtf.mu,linestyle = '--',label = '0.75mu')
    plt.legend()
    plt.xticks(rotation=90)
    plt.tight_layout()

    # plot2

    fig = plt.figure(figsize=(15, 3))
    plt.plot(x, dtf.evi_t1_t)
    plt.axhline(y=c, color='y', linestyle='--', label='c')
    # plt.axhline(y = c2, color = 'r', linestyle = '--', label='c2')
    plt.legend()
    plt.ylim([-1, 1])
    plt.xlabel('year_week_ts')
    plt.ylabel('EVI_t-1,t')
    plt.xticks(rotation=90)
    plt.tight_layout()

    # plot3

    fig = plt.figure(figsize=(15, 4))
    sns.scatterplot(x="year_week_ts", y="tseries", data=dtf, hue="ind",
                    palette="RdYlBu_r")  # gist_gray, flare,RdYlBu_r
    plt.xticks(rotation=90)

    # plot4

    fig = plt.figure(figsize=(15, 4))
    sns.scatterplot(x="year_week_ts", y="tseries", data=dtf, hue="ind2",
                    palette="RdYlBu_r")  # gist_gray, flare,RdYlBu_r
    plt.xticks(rotation=90)
