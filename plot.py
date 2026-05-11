import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
df = pd.read_csv('meta_plot/res_con.csv')
r = pd.read_csv('meta_plot/r_con.csv')

def plot_rank_across_rho(r):
    r[['diff_AY','diff_AMY','diff_AZY','diff_AZMY']] = np.where(r[['AY_ne','AMY_ne','AZY_ne','AZMY_ne']].to_numpy()-r[['AY_re','AMY_re','AZY_re','AZMY_re']].to_numpy()!=0,1,0)
    r['n_diff'] = r[['diff_AY','diff_AMY','diff_AZY','diff_AZMY']].sum(axis=1)
    unique_rho = r['rho'].unique()
    unique_a1 = r['a1'].unique()
    fig = plt.figure(figsize=(16, 12))
    for i, rho_value in enumerate(unique_rho):
        ax = fig.add_subplot(2, 2, i + 1, projection='3d')
        
        r_filtered = r[r['rho'] == rho_value]
        
        sc = ax.scatter(
            r_filtered['a1'], r_filtered['g6'], r_filtered['beta2'], 
            c=r_filtered['n_diff'], cmap='inferno', alpha=1.0
        )

        for a1_value in unique_a1:
            X, Y = np.meshgrid([a1_value], np.linspace(r_filtered['g6'].min(), r_filtered['g6'].max(), 100))
            Z = np.linspace(r_filtered['beta2'].min(), r_filtered['beta2'].max(), 100)
            Z = np.tile(Z, (len(Y), 1))
            ax.plot_surface(X, Y, Z, color='gray', alpha=0.3, edgecolor='none')
        
        cbar = plt.colorbar(sc, ax=ax, pad=0.1)
        cbar.set_label('n_diff')
        
        ax.set_title(f'rho = {rho_value}')
        ax.set_xlabel('a1')
        ax.set_ylabel('g6')
        ax.set_zlabel('beta2')
        ax.grid(True)
    plt.tight_layout()
    plt.savefig(f'meta_plot/rank_con.png', bbox_inches='tight')
    plt.show()
plot_rank_across_rho(r)

dfs = df
dfs['AZY_sign'] = np.where(dfs['AZY_ne'].to_numpy()*dfs['AZY_re'].to_numpy()>0,0,1)
dfs['AZMY_sign'] = np.where(dfs['AZMY_ne'].to_numpy()*dfs['AZMY_re'].to_numpy()>0,0,1)
print(dfs)
print(np.unique(dfs['AZY_sign']))
print(np.unique(dfs['AZMY_sign']))

def plot_sign_across_rho(data, sign_column):
    unique_rho = data['rho'].unique()
    unique_a1 = data['a1'].unique()
    fig = plt.figure(figsize=(16, 12))
    for i, rho_value in enumerate(unique_rho):
        ax = fig.add_subplot(2, 2, i + 1, projection='3d')
        
        data_filtered = data[data['rho'] == rho_value]
        
        sc = ax.scatter(
            data_filtered['a1'], data_filtered['g6'], data_filtered['beta2'], 
            c=data_filtered[sign_column], cmap='coolwarm', alpha=0.8
        )
        
        for a1_value in unique_a1:
            X, Y = np.meshgrid([a1_value], np.linspace(data_filtered['g6'].min(), data_filtered['g6'].max(), 10))
            Z = np.linspace(data_filtered['beta2'].min(), data_filtered['beta2'].max(), 10)
            Z = np.tile(Z, (len(Y), 1))
            ax.plot_surface(X, Y, Z, color='gray', alpha=0.3, edgecolor='none')
        
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label=f"{sign_column} same"),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label=f"{sign_column} opposite")
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        sign_counts = data_filtered[sign_column].value_counts()
        for sign_value, count in sign_counts.items():
            if sign_value==1:
                sign_text = f'opposite'
            else:
                sign_text = f'same'
            ax.text2D(0.05, 0.95 - 0.05 * sign_value, f"{sign_text}: {count}", transform=ax.transAxes)
        
        ax.set_title(f'rho = {rho_value}')
        ax.set_xlabel('a1')
        ax.set_ylabel('g6')
        ax.set_zlabel('beta2')
        ax.grid(True)
    plt.tight_layout()
    plt.savefig(f'meta_plot/{sign_column}_con.png', bbox_inches='tight')
    plt.show()
    
plot_sign_across_rho(dfs, 'AZY_sign')
plot_sign_across_rho(dfs, 'AZMY_sign')