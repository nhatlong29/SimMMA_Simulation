import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as plc
import pandas as pd
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm
from skimage.measure import marching_cubes
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_rank_across_rho(z, m, y, typ, sz):
    data = pd.read_csv(f'res_Z{z}_M{m}_Y{y}{typ}.csv')
    data['tol_re'] = data[['AY','AMY','AZY_re','AZMY_re','Z_e']].sum(axis=1)
    data['per_Ze'] = np.abs(data['Z_e'])/np.abs(data['tol_re'])
    data['mag_re'] = np.sum(np.abs(data[['AY','AMY','AZY_re','AZMY_re','Z_e']]),axis=1)
    data['per_magZe'] = np.abs(data['Z_e'])/np.abs(data['mag_re'])
    
    r = pd.read_csv(f'r_Z{z}_M{m}_Y{y}{typ}.csv')
    r = pd.merge(r,data[['per_Ze','per_magZe']],left_index=True,right_index=True)
    r['Z_e'] = np.abs(r['Z_e'])
    if isinstance(sz,str):
        r[sz] = r[sz].fillna(0)
        mask = (r[sz]==0)
        r.loc[~mask,'s'] = 6 + ((r.loc[~mask,sz] - r.loc[~mask,sz].min()) / (r.loc[~mask,sz].max() - r.loc[~mask,sz].min())) * (80 - 6)
        r.loc[mask,'s'] = 0
    else:
        r['s'] = 20
    print(r)
    r[['diff_AY','diff_AMY','diff_AZY','diff_AZMY']] = np.where(r[['AY_ne','AMY_ne','AZY_ne','AZMY_ne']].to_numpy()-r[['AY_re','AMY_re','AZY_re','AZMY_re']].to_numpy()!=0,1,0)
    r['n_diff'] = r[['diff_AY','diff_AMY','diff_AZY','diff_AZMY']].sum(axis=1)
    r.loc[r[['AY_ne','AMY_ne','AZY_ne','AZMY_ne','AY_re','AMY_re','AZY_re','AZMY_re']].isnull().values.any(axis=1),'n_diff'] = 1
    unique_rho = r['rho'].unique()
    unique_a1 = r['a1'].unique()
    colorlist = ['blue', 'gray', 'green', 'red', 'purple']
    fig = plt.figure(figsize=(16, 12))
    for i, rho_value in enumerate(unique_rho):
        ax = fig.add_subplot(1, 4, i + 1, projection='3d')
        
        r_filtered = r[r['rho'] == rho_value]
        custom_palette = ListedColormap(colorlist[:r_filtered['n_diff'].max()+1])
        sc = ax.scatter(
            r_filtered['a1'], r_filtered['g6'], r_filtered['beta2'], 
            c=r_filtered['n_diff'], cmap=custom_palette, norm=plc.Normalize(vmin=0, vmax=r_filtered['n_diff'].max(), clip=False), alpha=1.0,
            s=r_filtered['s']
        )
        for a1_value in unique_a1:
            X, Y = np.meshgrid([a1_value], np.linspace(r_filtered['g6'].min(), r_filtered['g6'].max(), 100))
            Z = np.linspace(r_filtered['beta2'].min(), r_filtered['beta2'].max(), 100)
            Z = np.tile(Z, (len(Y), 1))
            ax.plot_surface(X, Y, Z, color='silver', alpha=0.3, edgecolor='none')

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label=f"0"),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label=f"2"),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label=f"3"),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=10, label=f"4")
        ]
        #ax.legend(handles=legend_elements, loc='upper right')
        
        fig.legend(handles=legend_elements, loc='upper right')
        rank_counts = r_filtered['n_diff'].value_counts()
        for rank_value, count in rank_counts.items():
            rank_text = f'{rank_value} differences'
            ax.text2D(0.05, 0.95 - 0.05 * rank_value, f"{rank_text}: {count}", transform=ax.transAxes)
        
        xticks = ax.get_xticks()
        xticks_new = [tick for tick in xticks if tick in np.unique(r_filtered['a1'])]
        ax.set_xticks(xticks_new)
        yticks = ax.get_yticks()
        yticks_new = [tick for tick in yticks if tick in np.unique(r_filtered['g6'])]
        ax.set_yticks(yticks_new)
        zticks = ax.get_zticks()
        zticks_new = [tick for tick in zticks if tick in np.unique(r_filtered['beta2'])]
        ax.set_zticks(zticks_new)
        
        ax.set_title(f'rho = {rho_value}')
        ax.set_xlabel('alpha 1')
        ax.set_ylabel('gamma 6')
        ax.set_zlabel('beta 2')
        ax.grid(True)
        ax.view_init(elev=23, azim=-77)
    plt.tight_layout()
    #plt.savefig(f'rank_{sz}_Z{z}_M{m}_Y{y}{typ}.jpeg', bbox_inches='tight', dpi=900)
    plt.show()
#plot_rank_across_rho(z = 'con', m='con', y='con', typ='', sz = 'Z_e')
def plot_sign_across_rho(z, m, y, sign_column, typ, sz=None):
    data = pd.read_csv(f'res_Z{z}_M{m}_Y{y}{typ}.csv')
    
    data['tol_re'] = data[['AY','AMY','AZY_re','AZMY_re','Z_e']].sum(axis=1)
    data['per_Ze'] = np.abs(data['Z_e'])/np.abs(data['tol_re'])
    data['mag_re'] = np.sum(np.abs(data[['AY','AMY','AZY_re','AZMY_re','Z_e']]),axis=1)
    data['per_magZe'] = np.abs(data['Z_e'])/np.abs(data['mag_re'])
    
    data['AZY_sign'] = np.where(data['AZY_ne'].to_numpy()*data['AZY_re'].to_numpy()>0,0,1)
    mask = data[['AZY_ne', 'AZY_re']].isnull().values.any(axis=1)
    data.loc[mask, 'AZY_sign'] = 3
    data['AZMY_sign'] = np.where(data['AZMY_ne'].to_numpy()*data['AZMY_re'].to_numpy()>0,0,1)
    mask = data[['AZMY_ne', 'AZMY_re']].isnull().values.any(axis=1)
    data.loc[mask, 'AZMY_sign'] = 3
    
    data['Z_e'] = np.abs(data['Z_e'])
    if isinstance(sz,str):
        data[sz] = data[sz].fillna(0)
        mask = (data[sz]==0)
        data.loc[~mask,'s'] = 6 + ((data.loc[~mask,sz] - data.loc[~mask,sz].min()) / (data.loc[~mask,sz].max() - data.loc[~mask,sz].min())) * (60 - 6)
        data.loc[mask,'s'] = 0
    else:
        data['s'] = 20
    print(data)
    unique_rho = data['rho'].unique()
    unique_a1 = data['a1'].unique()
    fig = plt.figure(figsize=(16, 12))
    colorlist = ['blue', 'red', 'silver']
    for i, rho_value in enumerate(unique_rho):
        ax = fig.add_subplot(1, 4, i + 1, projection='3d')
        
        data_filtered = data[data['rho'] == rho_value]
        custom_palette = ListedColormap(colorlist[:data[sign_column].max()+1])
        sc = ax.scatter(
            data_filtered['a1'], data_filtered['g6'], data_filtered['beta2'], 
            c=data_filtered[sign_column], cmap=custom_palette, norm=plc.Normalize(vmin=0, vmax=data[sign_column].max(), clip=False), alpha=1.0,
            s=data_filtered['s']
        )
        
        for a1_value in unique_a1:
            X, Y = np.meshgrid([a1_value], np.linspace(data_filtered['g6'].min(), data_filtered['g6'].max(), 10))
            Z = np.linspace(data_filtered['beta2'].min(), data_filtered['beta2'].max(), 10)
            Z = np.tile(Z, (len(Y), 1))
            ax.plot_surface(X, Y, Z, color='silver', alpha=0.3, edgecolor='none')
        
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label=f"{sign_column} same"),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label=f"{sign_column} opposite")
        ]
        #ax.legend(handles=legend_elements, loc='upper right')
        
        fig.legend(handles=legend_elements, loc='upper right')
        sign_counts = data_filtered[sign_column].value_counts()
        for sign_value, count in sign_counts.items():
            if sign_value==1:
                ax.text2D(0.05, 0.95 - 0.05 * sign_value, f"Opposite: {count}", transform=ax.transAxes)
            elif sign_value==0:
                ax.text2D(0.05, 0.95 - 0.05 * sign_value, f"Same: {count}", transform=ax.transAxes)
            
        
        xticks = ax.get_xticks()
        xticks_new = [tick for tick in xticks if tick in np.unique(data_filtered['a1'])]
        ax.set_xticks(xticks_new)
        yticks = ax.get_yticks()
        yticks_new = [tick for tick in yticks if tick in np.unique(data_filtered['g6'])]
        ax.set_yticks(yticks_new)
        zticks = ax.get_zticks()
        zticks_new = [tick for tick in zticks if tick in np.unique(data_filtered['beta2'])]
        ax.set_zticks(zticks_new)      
        
        ax.set_title(f'rho = {rho_value}')
        ax.set_xlabel('alpha 1')
        ax.set_ylabel('gamma 6')
        ax.set_zlabel('beta 2')
        ax.grid(True)
        ax.view_init(elev=23, azim=-77)
    plt.tight_layout()
    #plt.savefig(f'up/{sign_column}_{sz}_Z{z}_M{m}_Y{y}{typ}.jpeg', bbox_inches='tight', dpi=900)
    plt.show()

def plot_combined(z, m, y, sz, ver, typ=''):
    data = pd.read_csv(f'res_Z{z}_M{m}_Y{y}{typ}.csv')
    data['tol_re'] = data[[f'P1_{ver}',f'P2_{ver}',f'P3_{ver}',f'P4_{ver}',f'Ze_{ver}']].sum(axis=1)
    data['per_Ze'] = np.abs(data[f'Ze_{ver}'])/np.abs(data['tol_re'])
    data['mag_re'] = np.sum(np.abs(data[[f'P1_{ver}',f'P2_{ver}',f'P3_{ver}',f'P4_{ver}',f'Ze_{ver}']]),axis=1)
    data['per_magZe'] = np.abs(data[f'Ze_{ver}'])/np.abs(data['mag_re'])
    #if typ=='':
    r = pd.read_csv(f'r_Z{z}_M{m}_Y{y}{typ}.csv')
    r = pd.merge(r,data[['per_Ze','per_magZe']],left_index=True,right_index=True)
    r[f'Ze_{ver}'] = np.abs(r[f'Ze_{ver}'])
    #else:
    #    r = pd.read_csv(f'r_Z{z}_M{m}_Y{y}{typ}.csv')
    #    r = pd.merge(r,data[['per_Ze','per_magZe',f'diff_P2_{ver}']],left_index=True,right_index=True)
    #    r[f'Ze_{ver}'] = np.abs(r[f'Ze_{ver}'])
    #    r[f'diff_P2_{ver}'] = np.abs(r[f'diff_P2_{ver}'])
    if isinstance(sz,str):
        r[f'{sz}_{ver}'] = r[f'{sz}_{ver}'].fillna(0)
        mask = (r[f'{sz}_{ver}']==0)
        r.loc[~mask,'s'] = 6 + ((r.loc[~mask,f'{sz}_{ver}'] - r.loc[~mask,f'{sz}_{ver}'].min()) / (r.loc[~mask,f'{sz}_{ver}'].max() - r.loc[~mask,f'{sz}_{ver}'].min())) * (90 - 6)
        r.loc[mask,'s'] = 3
    else:
        r['s'] = 20
    print(r)
    r[['rdiff_P1','rdiff_P2','rdiff_P3','rdiff_P4']] = np.where(r[['P1','P2','P3','P4']].to_numpy()-r[[f'P1_{ver}',f'P2_{ver}',f'P3_{ver}',f'P4_{ver}']].to_numpy()!=0,1,0)
    r[['rdiff_val_P1','rdiff_val_P2','rdiff_val_P3','rdiff_val_P4']] = np.where(r[['rdiff_P1','rdiff_P2','rdiff_P3','rdiff_P4']].to_numpy()!=0,
                                                                            data[['P1','P2','P3','P4']].to_numpy()-data[[f'P1_{ver}',f'P2_{ver}',f'P3_{ver}',f'P4_{ver}']].to_numpy(),0)
    r['n_rdiff'] = r[['rdiff_P1','rdiff_P2','rdiff_P3','rdiff_P4']].sum(axis=1)
    r['val_rdiff'] = r[['rdiff_val_P1','rdiff_val_P2','rdiff_val_P3','rdiff_val_P4']].abs().sum(axis=1)
    r.loc[r[['P1','P2','P3','P4',f'P1_{ver}',f'P2_{ver}',f'P3_{ver}',f'P4_{ver}']].isnull().values.any(axis=1),'n_rdiff'] = 1
    unique_rho = r['rho'].unique()
    unique_a1 = r['a1'].unique()
    colorlist = ['blue', 'gray', 'green', 'red', 'purple']
    fig, axes = plt.subplots(3, 4, figsize=(24, 18), subplot_kw={'projection': '3d'})
    # 1st row: ranking
    for i, rho_value in enumerate(unique_rho):
        ax = axes[0,i]
        r_filtered = r[r['rho'] == rho_value]
        custom_palette = ListedColormap(colorlist[:r_filtered['n_rdiff'].max()+1])
        sc = ax.scatter(
            r_filtered['a1'], r_filtered['g6'], r_filtered['beta2'], 
            c=r_filtered['n_rdiff'], cmap=custom_palette, norm=plc.Normalize(vmin=0, vmax=r_filtered['n_rdiff'].max(), clip=False), alpha=1.0,
            s=r_filtered['s']
        )
        for a1_value in unique_a1:
            X, Y = np.meshgrid([a1_value], np.linspace(r_filtered['g6'].min(), r_filtered['g6'].max(), 100))
            Z = np.linspace(r_filtered['beta2'].min(), r_filtered['beta2'].max(), 100)
            Z = np.tile(Z, (len(Y), 1))
            ax.plot_surface(X, Y, Z, color='silver', alpha=0.3, edgecolor='none')

        xticks = ax.get_xticks()
        xticks_new = [tick for tick in xticks if tick in np.unique(r_filtered['a1'])]
        ax.set_xticks(xticks_new)
        yticks = ax.get_yticks()
        yticks_new = [tick for tick in yticks if tick in np.unique(r_filtered['g6'])]
        ax.set_yticks(yticks_new)
        zticks = ax.get_zticks()
        zticks_new = [tick for tick in zticks if tick in np.unique(r_filtered['beta2'])]
        ax.set_zticks(zticks_new)
        ax.set_xlabel(r'$\alpha_1$',fontsize=10)
        ax.set_ylabel(r'$\gamma_6$',fontsize=10)
        ax.set_zlabel(r'$\beta_2$',fontsize=10)
        ax.grid(True)
        ax.view_init(elev=23, azim=-77)
    
    data[f'P2_sign'] = np.where(data['P2'].to_numpy()*data[f'P2_{ver}'].to_numpy()>0,0,1)
    mask = data[['P2', f'P2_{ver}']].isnull().values.any(axis=1)
    data.loc[mask, f'P2_sign'] = 3
    data[f'P3_sign'] = np.where(data['P3'].to_numpy()*data[f'P3_{ver}'].to_numpy()>0,0,1)
    mask = data[['P3', f'P3_{ver}']].isnull().values.any(axis=1)
    data.loc[mask, f'P3_sign'] = 3
    data[f'Ze_{ver}'] = np.abs(data[f'Ze_{ver}'])
    #data[f'diff_P2_{ver}'] = np.abs(data[f'diff_P2_{ver}'])
    if isinstance(sz,str):
        data[f'{sz}_{ver}'] = data[f'{sz}_{ver}'].fillna(0)
        mask = (data[f'{sz}_{ver}']==0)
        data.loc[~mask,'s'] = 6 + ((data.loc[~mask,f'{sz}_{ver}'] - data.loc[~mask,f'{sz}_{ver}'].min()) / (data.loc[~mask,f'{sz}_{ver}'].max() - data.loc[~mask,f'{sz}_{ver}'].min())) * (90 - 6)
        data.loc[mask,'s'] = 3
    else:
        data['s'] = 20
    print(data)
    unique_rho = data['rho'].unique()
    unique_a1 = data['a1'].unique()
    colorlist = ['blue', 'red', 'gray']
    for j, sign in enumerate(['P2', 'P3']):    
        for i, rho_value in enumerate(unique_rho):
            sign_column = f'{sign}_sign'
            ax = axes[j+1, i]
            data_filtered = data[data['rho'] == rho_value]
            custom_palette = ListedColormap(colorlist[:data[sign_column].max()+1])
            sc = ax.scatter(
                data_filtered['a1'], data_filtered['g6'], data_filtered['beta2'], 
                c=data_filtered[sign_column], cmap=custom_palette, norm=plc.Normalize(vmin=0, vmax=data[sign_column].max(), clip=False), alpha=1.0,
                s=data_filtered['s']
            )
            for a1_value in unique_a1:
                X, Y = np.meshgrid([a1_value], np.linspace(data_filtered['g6'].min(), data_filtered['g6'].max(), 10))
                Z = np.linspace(data_filtered['beta2'].min(), data_filtered['beta2'].max(), 10)
                Z = np.tile(Z, (len(Y), 1))
                ax.plot_surface(X, Y, Z, color='silver', alpha=0.3, edgecolor='none')
                        
            xticks = ax.get_xticks()
            xticks_new = [tick for tick in xticks if tick in np.unique(data_filtered['a1'])]
            ax.set_xticks(xticks_new)
            yticks = ax.get_yticks()
            yticks_new = [tick for tick in yticks if tick in np.unique(data_filtered['g6'])]
            ax.set_yticks(yticks_new)
            zticks = ax.get_zticks()
            zticks_new = [tick for tick in zticks if tick in np.unique(data_filtered['beta2'])]
            ax.set_zticks(zticks_new)  
            ax.set_xlabel(r'$\alpha_1$',fontsize=12)
            ax.set_ylabel(r'$\gamma_6$',fontsize=12)
            ax.set_zlabel(r'$\beta_2$',fontsize=12)
            ax.grid(True)
            ax.view_init(elev=23, azim=-77)
    legend2 = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label=f"Same sign"),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label=f"Opposite sign")
        ]
    legend1 = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label=f"0"),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label=f"2"),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label=f"3"),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=10, label=f"4")
        ]
    axes[0,0].legend(handles=legend1, loc='upper left', fontsize=12)
    axes[1,0].legend(handles=legend2, loc='upper left', fontsize=12)
    axes[2,0].legend(handles=legend2, loc='upper left', fontsize=12)
    # Set column titles (rho values) on top
    for i, rho_value in enumerate(unique_rho):
        axes[0, i].set_title(f"rho = {rho_value}", fontsize=16, fontweight='bold')
    row_titles = ["Rank", "A → Z → Y sign", "A → Z → M → Y sign"]
    for row, title in enumerate(row_titles):
        fig.text(0.18, 0.8 - row * 0.3, title, va='center', ha='right', rotation='vertical', fontsize=14, fontweight='bold')
    for ax in axes.flat:
        ax.tick_params(axis='both', labelsize=12)
    plt.tight_layout(rect=[0.1, 0, 1, 1], pad=1.0, w_pad=-1, h_pad=0)
    plt.subplots_adjust(left=0.12, right=0.98, top=0.92, bottom=0.08, wspace=-0.3, hspace=0.15)
    plt.savefig(f'up/{sz}_Z{z}_M{m}_Y{y}{typ}_ver{ver}.jpeg', bbox_inches='tight', dpi=900)
    #plt.show()
    # export data
    df_out = pd.merge(left=data, right=r[['a1','g6','beta2','rho','n_rdiff','val_rdiff']], how='outer', on=['a1','g6','beta2','rho'])
    df_out['Z'] = z
    df_out['M'] = m
    df_out['Y'] = y
    df_out['ver'] = ver
    return df_out

#df = plot_combined(z='con', m='con', y='bi', typ='upbound', sz='diff_P2', ver='1.0')
# Usage

"""sz = "Ze"
df_all = pd.DataFrame()
for z in ['bi','con']:
    for m in ['bi','con']:
        for y in ['bi','con']:
            for ver in ['1.0']: #,'2.0','3.0'
                df = plot_combined(z=z, m=m, y=y, sz=sz, ver=ver)
                #df_all = pd.concat([df_all,df], axis=0, ignore_index=True)"""

def kount(data, sub, id, col):
    if not sub:
        subset = data
    else:
        subset = data[np.logical_and.reduce([
        data[k].isin(v) if isinstance(v, (list, tuple, set, np.ndarray)) else data[k] == v
        for k, v in sub.items()])]
    df_melt = subset.melt(id_vars=id, value_vars=col,
                  var_name='col', value_name='val')
    df_melt['col_val'] = df_melt['col'] + '_' + df_melt['val'].astype(str)
    group_col_name = '_'.join(id)
    df_melt[group_col_name] = df_melt[id].astype(str).agg('_'.join, axis=1)
    matrix = pd.crosstab(df_melt[group_col_name], df_melt['col_val'])
    matrix = matrix.sort_values(by=group_col_name)
    return matrix


def count_rank_and_sign_across_rho(typ, ver): # need to be updated.....
    for mod in ['rank','sign']:
        df_all = pd.DataFrame()
        if mod == 'sign':
            for z in ['bi','con']:
                for m in ['bi','con']:
                    for y in ['bi','con']:                       
                        df = pd.read_csv(f'res_Z{z}_M{m}_Y{y}{typ}.csv')
                        df['Z'] = z
                        df['M'] = m
                        df['Y'] = y
                        df_all = pd.concat([df_all,df], axis=0, ignore_index=True)
            df_all[f'P2_sign'] = np.where(df_all['P2'].to_numpy()*df_all[f'P2_{ver}'].to_numpy()>0,0,1)
            mask = df_all[['P2', f'P2_{ver}']].isnull().values.any(axis=1)
            df_all.loc[mask, f'P2_sign'] = 3
            df_all[f'P3_sign'] = np.where(df_all['P3'].to_numpy()*df_all[f'P3_{ver}'].to_numpy()>0,0,1)
            mask = df_all[['P3', f'P3_{ver}']].isnull().values.any(axis=1)
            df_all.loc[mask, f'P3_sign'] = 3
            for sub in [{},
                        {'Z':'con','M':'con','Y':'con'},
                        {'Z':'con','M':'con','Y':'bi'},
                        {'Z':'con','M':'bi','Y':'con'},
                        {'Z':'con','M':'bi','Y':'bi'},
                        {'Z':'bi','M':'con','Y':'con'},
                        {'Z':'bi','M':'con','Y':'bi'},
                        {'Z':'bi','M':'bi','Y':'con'},
                        {'Z':'bi','M':'bi','Y':'bi'}]:
                if not sub:
                    id = ['Z','M','Y']
                else:
                    id = ['a1']
                test = kount(df_all, sub=sub, id=id, col=['P2_sign','P3_sign'])
                with open(f'up/sign.txt', 'a') as f:
                    f.write(str(sub) + '\n')
                    f.write(test.to_string() + '\n')
                    f.write('\n')
        elif mod == 'rank':
            for z in ['bi','con']:
                for m in ['bi','con']:
                    for y in ['bi','con']:                       
                        data = pd.read_csv(f'r_Z{z}_M{m}_Y{y}{typ}.csv')
                        data['Z'] = z
                        data['M'] = m
                        data['Y'] = y
                        df_all = pd.concat([df_all,data], axis=0, ignore_index=True)
                        
            df_all[['diff_P1','diff_P2','diff_P3','diff_P4']] = np.where(df_all[['P1','P2','P3','P4']].to_numpy()-df_all[[f'P1_{ver}',f'P2_{ver}',f'P3_{ver}',f'P4_{ver}']].to_numpy()!=0,1,0)
            df_all['n_diff'] = df_all[['diff_P1','diff_P2','diff_P3','diff_P4']].sum(axis=1)
            df_all.loc[df_all[['P1','P2','P3','P4',f'P1_{ver}',f'P2_{ver}',f'P3_{ver}',f'P4_{ver}']].isnull().values.any(axis=1),'n_diff'] = 1

            for sub in [{},
                        {'Z':'con','M':'con','Y':'con'},
                        {'Z':'con','M':'con','Y':'bi'},
                        {'Z':'con','M':'bi','Y':'con'},
                        {'Z':'con','M':'bi','Y':'bi'},
                        {'Z':'bi','M':'con','Y':'con'},
                        {'Z':'bi','M':'con','Y':'bi'},
                        {'Z':'bi','M':'bi','Y':'con'},
                        {'Z':'bi','M':'bi','Y':'bi'}]:
                if not sub:
                    id = ['Z','M','Y']
                else:
                    id = ['a1']
                test = kount(df_all, sub=sub, id=id, col=['n_diff'])
                with open(f'up/rank.txt', 'a') as f:
                    f.write(str(sub) + '\n')
                    f.write(test.to_string() + '\n')
                    f.write('\n')
#count_rank_and_sign_across_rho(typ='', ver="1.0")
def count_bound():
    df_all = pd.DataFrame()
    for z in ['bi','con']:
        for m in ['bi','con']:
            for y in ['bi']:                       
                df = pd.read_csv(f'res_Z{z}_M{m}_Y{y}upbound.csv')
                df['Z'] = z
                df['M'] = m
                df['Y'] = y
                df_all = pd.concat([df_all,df], axis=0, ignore_index=True)
       
    df_all['abs_diff'] = np.abs(df_all['diff_P2_1.0'])
    df_all['upboundZ'] = df_all['uppb_P2_z_1.0']
    df_all['upbound1'] = df_all['uppb_P2_1_1.0']
    df_all['upboundbo'] = df_all['uppb_P2_bo_1.0']
    df_all['upboundbn'] = df_all['uppb_P2_bn_1.0']
    df_all['upboundopt'] = df_all[['upboundZ','upbound1','upboundbo','upboundbn']].min(axis=1)
    
    df_all['lwboundZ'] = df_all['diff_PZ_1.0'] - df_all['uppb_P3_z_1.0']
    df_all['lwbound1'] = df_all['diff_PZ_1.0'] - df_all['uppb_P3_1_1.0']
    df_all['lwboundbo'] = df_all['diff_PZ_1.0'] - df_all['uppb_P3_bo_1.0']
    df_all['lwboundbn'] = df_all['diff_PZ_1.0'] - df_all['uppb_P3_bn_1.0']
    df_all['lwboundopt'] = df_all['diff_PZ_1.0'] - df_all[['uppb_P3_z_1.0','uppb_P3_1_1.0','uppb_P3_bo_1.0','uppb_P3_bn_1.0']].min(axis=1)
    
    cols = ['abs_diff', 'upboundZ', 'upbound1', 'upboundbo', 'upboundbn', 'upboundopt', 'lwboundZ', 'lwbound1', 'lwboundbo', 'lwboundbn', 'lwboundopt']
    df_all['checkup'] = df_all[['upboundZ','upbound1','upboundbn']].idxmin(axis=1)
    df_all['checklw'] = df_all[['lwboundZ','lwbound1','lwboundbn']].idxmax(axis=1)
    for z in ['bi','con']:
        for m in ['bi','con']:
            print(len(df_all[(df_all['Z']==z)&(df_all['M']==m)]))
            print(f'z {z} m {m}: {df_all.loc[((df_all['Z']==z)&(df_all['M']==m)),'checkup'].value_counts()}')
            print(f'z {z} m {m}: {df_all.loc[((df_all['Z']==z)&(df_all['M']==m)),'checklw'].value_counts()}')
    df_all["group"] = ("Z " + df_all["Z"].map({"bi": "binary", "con": "continuous"})
                       + " M " + df_all["M"].map({"bi": "binary", "con": "continuous"}))
    df_all = df_all.reset_index(drop=True)
    df_all["id"] = df_all.index
    # reshape wide -> long
    df_long = df_all.melt(
        id_vars=["id", "group"],
        value_vars=cols,
        var_name="Variable",
        value_name="Value"
    )
    groups = df_long["group"].unique()
    fig, axes = plt.subplots(int((len(cols)-1)/2), len(groups), figsize=(24, 3*(len(cols)-1)/2), sharey=False)
    thr = df_all[cols[0]].median()
    for tar in range(int((len(cols)-1)/2)):
        for j, grp in enumerate(groups):
            ax = axes[tar, j]
            sub = df_long[df_long["group"] == grp]
            wide = sub.pivot(index="id", columns="Variable", values="Value")
            abs_vals = wide[cols[0]]
            for i, row in wide.iterrows():
                color = "red" if abs_vals.loc[i] < thr else "blue"
                ax.plot(
                    [cols[tar+int((len(cols)-1)/2+1)],cols[0], cols[tar+1]],
                    [row[cols[tar+int((len(cols)-1)/2+1)]], row[cols[0]], row[cols[tar+1]]],
                    color=color,
                    alpha=0.2,
                    linewidth=1
                )
            ax.set_ylim(-1.7, 1.7)
            ax.scatter(wide.index.map(lambda _: cols[0]),
                    wide[cols[0]], s=10)
            ax.scatter(wide.index.map(lambda _: cols[tar+1]),
                    wide[cols[tar+1]], s=10)
            ax.scatter(wide.index.map(lambda _: cols[tar+int((len(cols)-1)/2+1)]),
                    wide[cols[tar+int((len(cols)-1)/2+1)]], s=10)
            ax.set_title(grp)
        axes[tar,0].set_ylabel(f"{cols[0]} → {cols[tar+1]}")
    plt.tight_layout()
    #plt.savefig(f'up/bound.jpeg', bbox_inches='tight', dpi=900)
    plt.show()
#count_bound()

def final_bound():
    df_all = pd.DataFrame()
    for z in ['bi','con']:
        for m in ['bi','con']:
            for y in ['bi']:
                df = pd.read_csv(f'res_Z{z}_M{m}_Y{y}upbound_upfin.csv')
                df['Z'] = z
                df['M'] = m
                df['Y'] = y
                df_all = pd.concat([df_all,df], axis=0, ignore_index=True)
       
    df_all['Absolute difference'] = np.abs(df_all['diff_P2_1.0'])
    df_all['upboundZ'] = df_all['uppb_P2_z_1.0']
    df_all['upbound1'] = df_all['uppb_P2_1_1.0']
    df_all['upboundbo'] = df_all['uppb_P2_bo_1.0']
    df_all['upboundbn'] = df_all['uppb_P2_bn_1.0']
    df_all['Optimal upper bound'] = df_all[['upboundZ','upbound1','upboundbo','upboundbn']].min(axis=1)
        
    cols = ['Absolute difference', 'Optimal upper bound']
    df_all["group"] = (df_all["Z"].map({"bi": "Binary", "con": "Continuous"}) + " Z, " 
                       + df_all["M"].map({"bi": "binary", "con": "continuous"}) + " M" )
    df_all = df_all.reset_index(drop=True)
    df_all["id"] = df_all.index
    # reshape wide -> long
    df_long = df_all.melt(
        id_vars=["id", "group"],
        value_vars=cols,
        var_name="Variable",
        value_name="Value"
    )
    groups = df_long["group"].unique()
    fig, axes = plt.subplots(2, 2, figsize=(16, 8), sharey=False)
    thr = df_all[cols[0]].median()
    
    for j, grp in enumerate(groups):
        ax = axes[0 if j<=1 else 1, 0 if j%2==0 else 1]
        sub = df_long[df_long["group"] == grp]
        wide = sub.pivot(index="id", columns="Variable", values="Value")
        abs_vals = wide[cols[0]]
        for i, row in wide.iterrows():
            color = "red" if abs_vals.loc[i] < thr else "blue"
            ax.plot(
                [cols[0], cols[1]],
                [row[cols[0]], row[cols[1]]],
                color=color,
                alpha=0.2,
                linewidth=1
            )
        ax.set_ylim(-0.1, 1.0)
        ax.scatter(wide.index.map(lambda _: cols[0]),
                wide[cols[0]], s=10)
        ax.scatter(wide.index.map(lambda _: cols[1]),
                wide[cols[1]], s=10)
        ax.set_title(grp, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'up/optupbound_fin.jpeg', bbox_inches='tight', dpi=900)
    plt.show()
#final_bound()

def alpha_range(a0 = np.linspace(-5, 5, 100), a1 = np.linspace(-5, 5, 100), rho = np.linspace(-1, 1, 100)):    
    def f(a0, a1, rho, names : str):
        mu0 = np.exp(a0)/(1+np.exp(a0))
        mu1 = np.exp(a0+a1)/(1+np.exp(a0+a1))
        p0 = mu1 - rho*np.sqrt(mu0*mu1*(1-mu1)/(1-mu0))
        p1 = mu1 + rho*np.sqrt((1-mu0)*mu1*(1-mu1)/mu0)
        if names == 'p_z0_1':
            return p1
        elif names == 'p_z0_0':
            return p0
        else:
            mask = (p1 > 0) & (p1 < 1) & (p0 > 0) & (p0 < 1)
            result = np.full_like(p1, -1, dtype=float)
            result[mask] = p1[mask] 
            return result
    
    A0, A1, RHO = np.meshgrid(a0, a1, rho, indexing='ij')
    G = f(A0, A1, RHO, 'all')
    verts, faces, normals, values = marching_cubes(G, level=0, spacing=(a0[1]-a0[0], a1[1]-a1[0], rho[1]-rho[0]))
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_trisurf(verts[:, 0] + a0[0], verts[:, 1] + a1[0], verts[:, 2] + rho[0], triangles=faces, color='gray', alpha=0.5)

    ax.set_xlabel('a0')
    ax.set_ylabel('a1')
    ax.set_zlabel('rho')
    #ax.set_title('Available region')
    ax.view_init(elev=50, azim=-60)
    #plt.savefig(f'up/a0_a1_ ranging.jpeg', bbox_inches='tight', dpi=900)
    plt.show()

def diff_p2_range_ver0(s1_1=np.linspace(0, 1, 100),
                       s1_prime_1=np.linspace(0, 1, 100),
                       s2_prime_1=np.linspace(0, 1, 100),
                       s3_1=np.linspace(0, 1, 100)):
    def f(a, b, c, d):
        part0 = np.abs(np.maximum(-1,a-d-1)-(b-c))
        part1 = np.abs(np.minimum(1,a-d+1)-(b-c))
        return np.maximum(part0,part1)
    
    S1, S2, S3, S4 = np.meshgrid(s1_1, s1_prime_1, s2_prime_1, s3_1, indexing='ij')
    G = f(S1, S2, S3, S4)
    
    slices = [0, 0.5, 1]
    idx = [np.argmin(np.abs(s3_1 - s)) for s in slices]
    
    fig = plt.figure(figsize=(18, 6))
    
    # normalize for coloring
    norm = plt.Normalize(G.min(), G.max())
    cmap = 'viridis'
    for i, (s_val, k) in enumerate(zip(slices, idx)):

        ax = fig.add_subplot(1, 3, i+1, projection='3d')

        G_slice = G[:, :, :, k]

        sc = ax.scatter(
            S1[:, :, :, k].flatten(),
            S2[:, :, :, k].flatten(),
            S3[:, :, :, k].flatten(),
            c=G_slice.flatten(),
            cmap=cmap,
            norm=norm,
            s=3,
            alpha=0.4
        )

        ax.set_title(f"s3_1 = {s_val}")
        ax.set_xlabel("s1_1")
        ax.set_ylabel("s1'_1")
        ax.set_zlabel("s2'_1")
        ax.view_init(elev=10, azim=-30)
        
    cbar = fig.colorbar(sc, ax=fig.axes, shrink=0.6, pad=0.04)
    cbar.set_label("bound value")
    plt.tight_layout()
    plt.show()
#diff_p2_range_ver0()

def diff_p2_range_ver1(typ, s1_1=np.linspace(0, 1, 100),
                       s1_prime_1=np.linspace(0, 1, 100),
                       s2_prime_1=np.linspace(0, 1, 100)):
    def f(a, b, c):
        part0 = np.maximum.reduce([a, 1 - b, c])
        part1 = np.minimum.reduce([a, 1 - b, c])
        part2 = np.abs(b - (1 - c))
        return 2 - (part0 - part1) - part2
    def f1(a, b, c):
        part0 = np.abs((a-1) - (b-c))
        part1 = np.abs(a - (b-c))
        return np.maximum(part0,part1)
    def f2(a, b, c):
        part0 = np.abs(np.maximum(-1,a-1) - (b-c))
        part1 = np.abs(np.minimum(1,a+1) - (b-c))
        return np.maximum(part0,part1)
    
    S1, S2, S3 = np.meshgrid(s1_1, s1_prime_1, s2_prime_1, indexing='ij')
    if typ == 0:
        G = f(S1, S2, S3)
    elif typ == 1:
        G = f1(S1, S2, S3)
    else:
        s1_1 = np.linspace(-1,1,100)
        S1, S2, S3 = np.meshgrid(s1_1, s1_prime_1, s2_prime_1, indexing='ij')
        G = f2(S1, S2, S3)
        
    fig = plt.figure(figsize=(18, 6))
    
    # normalize for coloring
    norm = plt.Normalize(0, 2)
    
    ax = fig.add_subplot(131, projection='3d')
    sc = ax.scatter(
            S1.flatten(),
            S2.flatten(),
            S3.flatten(),
            c=G.flatten(),
            norm=norm,
            cmap="viridis",
            alpha=0.2,
            s=5
        )
    
    ax.set_title('Full', pad=2)    
    
    ax1 = fig.add_subplot(132, projection='3d')
    mask_low = G < 1.0
    sc1 = ax1.scatter(
        S1[mask_low],
        S2[mask_low],
        S3[mask_low],
        c=G[mask_low],
        norm=norm,
        cmap="viridis",
        s=5,
        alpha=0.2
    )
    ax1.set_title("Bound < 1.0", pad=2)
    
    ax2 = fig.add_subplot(133, projection='3d')
    mask_high = G >= 1.0
    sc2 = ax2.scatter(
        S1[mask_high],
        S2[mask_high],
        S3[mask_high],
        c=G[mask_high],
        norm=norm,
        cmap="viridis",
        s=5,
        alpha=0.2
    )
    ax2.set_title("Bound ≥ 1.0", pad=2)

    cbar = fig.colorbar(sc, ax=[ax, ax1, ax2], shrink=0.6,  fraction=0.03, pad=0.04)
    cbar.set_label("Bound value")

    for ax in [ax, ax1, ax2]:
        if typ==0 or typ ==1:
            ax.set_xlabel("s1_1")
        else:
            ax.set_xlabel("psi_z")
        ax.set_ylabel("s1'_1")
        ax.set_zlabel("s2'_1")
        ax.view_init(elev=10, azim=-30)
        
    plt.subplots_adjust(wspace=0.05, right=0.85)
    plt.savefig(f'up/bound{typ}.jpeg', bbox_inches='tight', dpi=900)
    #plt.show()
#diff_p2_range_ver1(typ=2)
def bound_across_params_3d(path, axx='beta2', axy='a1', p='g2'):
    df_all = pd.DataFrame()
    for z in ['bi','con']:
        for m in ['bi','con']:
            for y in ['bi']:
                df_list = [pd.read_csv(f'res_Z{z}_M{m}_Y{y}upbound_newset.csv')] #pd.read_csv(f'res_Z{z}_M{m}_Y{y}upbound.csv'),
                for df in df_list:
                    df['Z'] = z
                    df['M'] = m
                    df['Y'] = y
                    df_all = pd.concat([df_all,df], axis=0, ignore_index=True)
    if path=='P3':                
        df_all[f'diff_P3_1.0'] = df_all['Ys2'] - df_all['Ys3'] - df_all['Ys2.dprime'] + df_all['Ys3.dprime']
    df_all['Absolute difference'] = np.abs(df_all[f'diff_{path}_1.0'])
    df_all['upboundZ'] = df_all[f'uppb_{path}_z_1.0']
    df_all['upbound1'] = df_all[f'uppb_{path}_1_1.0']
    df_all['upboundbo'] = df_all[f'uppb_{path}_bo_1.0']
    df_all['upboundbn'] = df_all[f'uppb_{path}_bn_1.0']
    df_all['Optimal upper bound'] = df_all[['upboundZ','upbound1','upboundbo','upboundbn']].min(axis=1)
    df_all['diff_bound_true'] = df_all['Optimal upper bound'] - df_all['Absolute difference']
    
    z_m_order = [('con', 'con'), ('con', 'bi'), ('bi', 'con'), ('bi', 'bi')]
    rho_values = sorted(df_all['rho'].unique())
    p_values = sorted(df_all[p].unique())
    cmap = plt.cm.get_cmap('viridis', len(p_values))
    
    fig, axes = plt.subplots(4, 4, figsize=(30, 24), subplot_kw={'projection': '3d'})
    for row, (z, m) in enumerate(z_m_order):
        for col, rho in enumerate(rho_values):
            ax = axes[row, col]
            df_filtered = df_all[(df_all['Z'] == z) & (df_all['M'] == m) & (df_all['rho'] == rho)]
            df_filtered = df_filtered.dropna(subset=[p, axx, axy, 'diff_bound_true'])
            
            if df_filtered.empty:
                ax.text2D(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center')
                #ax.set_title(f'Z={z}, M={m}, rho={rho}')
                continue
            
            
            for k, val in enumerate(p_values):
                plane = df_filtered[df_filtered[p] == val]
                if plane.empty:
                    continue
                pivot = plane.pivot(index=axy, columns=axx, values='diff_bound_true')
                if pivot.empty:
                    continue
                X, Y = np.meshgrid(pivot.columns.values, pivot.index.values)
                Z = pivot.values
                ax.plot_surface(X, Y, Z, color=cmap(k), alpha=0.45, edgecolor='none', rstride=1, cstride=1)
            
            xticks = ax.get_xticks()
            xticks_new = [tick for tick in xticks if tick in np.unique(df_filtered[axx])]
            ax.set_xticks(xticks_new)
            yticks = ax.get_yticks()
            yticks_new = [tick for tick in yticks if tick in np.unique(df_filtered[axy])]
            ax.set_yticks(yticks_new)
            zticks = ax.get_zticks()
            zticks_new = [round(tick,2) for tick in zticks]
            ax.set_zticks(zticks_new)
            ax.set_xlabel(r'$\beta_2$', fontsize=12)
            ax.set_ylabel(r'$\alpha_1$', fontsize=12)
            ax.set_zlabel(f'Difference', fontsize=12)
            #ax.set_title(f'Z={z}, M={m}, rho={rho}', fontsize=11)
            ax.view_init(elev=25, azim=-25)
    
    legend_elements = [Line2D([0], [0], marker='s', color='w', markerfacecolor=cmap(i), markersize=10, label=f'$\gamma_6$ = {p}') for i, p in enumerate(p_values)]
    # Set column titles (rho values) on top
    for i, rho_value in enumerate(np.unique(df_all['rho'])):
        axes[0, i].set_title(f"rho = {rho_value}", fontsize=16, fontweight='bold')
    row_titles = ["Z continuous, M continuous", "Z continuous, M binary", "Z binary, M continuous", "Z binary, M binary"]
    for row, title in enumerate(row_titles):
        fig.text(0.16, 0.83 - row * 0.22, title, va='center', ha='right', rotation='vertical', fontsize=14, fontweight='bold')
    for ax in axes.flat:
        ax.tick_params(axis='both', labelsize=12)
    fig.legend(handles=legend_elements, loc='upper right',  bbox_to_anchor=(0.95, 0.95),  fontsize=12)
    plt.tight_layout(rect=[0.1, 0, 1, 1], pad=1.0, w_pad=-1, h_pad=0)
    plt.subplots_adjust(left=0.15, right=0.93, top=0.95, bottom=0.05, wspace=0, hspace=0)
    plt.savefig(f'up/bound{path}_new.pdf', bbox_inches='tight')
    #plt.show()
#bound_across_params_3d(path="P2", axx='beta2', axy='a1', p='g2')

def bound_across_params_2d(conditions = {"rho":0.2, "g6":-1.5, "beta2":1.5}, free = 'a1'):
    df_all = pd.DataFrame()
    """for z in ['bi','con']:
        for m in ['bi','con']:
            for y in ['bi']:
                df_list = [pd.read_csv(f'res_Z{z}_M{m}_Y{y}upbound_newset.csv')] #pd.read_csv(f'res_Z{z}_M{m}_Y{y}upbound.csv'),
                for df in df_list:
                    df['Z'] = z
                    df['M'] = m
                    df['Y'] = y
                    df_all = pd.concat([df_all,df], axis=0, ignore_index=True)"""
    
    df_all = pd.read_csv(f'res_test.csv')
    df_all['Z'] = df_all['mode.z']
    df_all['M'] = df_all['mode.m']
    df_all['Y'] = df_all['mode.y']
    df_all = df_all[df_all['Y']=='bi']
    mask = pd.Series(True, index=df_all.index)
    for col, value in conditions.items():
        mask &= (df_all[col] == value)
    df_all = df_all[mask]
    # 3 bounds
    df_all[f'abs_diff_P2'] = np.abs(df_all[f'diff_P2_1.0'])
    df_all[f'abs_diff_P3'] = np.abs(df_all['Ys2'] - df_all['Ys3'] - df_all['Ys2.dprime'] + df_all['Ys3.dprime'])
    for p in ['P2','P3']:        
        df_all[f'boundopt_{p}'] = df_all[[f'uppb_{p}_z_1.0',f'uppb_{p}_1_1.0',f'uppb_{p}_bo_1.0',f'uppb_{p}_bn_1.0']].min(axis=1)   
    df_all['abs_diff_P2P3'] = df_all['abs_diff_P2'] + df_all['abs_diff_P3']
    df_all['boundopt_P2P3'] = np.abs(df_all['Ys1'] - df_all['Ys3'] - df_all['Ys1.prime'] + df_all['Ys3.dprime'])
    df_all['upboundopt_P2P3'] = df_all['boundopt_P2'] + df_all['boundopt_P3'] 
    z_m_order = [('con', 'con'), ('con', 'bi'), ('bi', 'con'), ('bi', 'bi')]
    # plot    
    if len(conditions)==3:
        fig, axes = plt.subplots(4, 3, figsize=(16, 12), sharey=False)
        for row, (z, m) in enumerate(z_m_order):
            for col, p in enumerate(['P2','P3','P2P3']):
                ax = axes[row, col]
                df_filtered = df_all[(df_all['Z'] == z) & (df_all['M'] == m)]
                df_filtered = df_filtered.sort_values(free)
                
                ax.plot(df_filtered[free], df_filtered[f'abs_diff_{p}'] , color='black', linewidth=1)
                ax.plot(df_filtered[free], df_filtered[f'boundopt_{p}'] , color='black', linestyle='--', linewidth=1)
                
                ax.set_ylim(-0.1, 1.0)
                if free=='a1' or free=='alpha_1':
                    ax.set_xlabel(r'$\alpha_1$')
                elif free=='g6' or free=='gamma_6':
                    ax.set_xlabel(r'$\gamma_6$')
                elif free=='g2' or free=='gamma_2':
                    ax.set_xlabel(r'$\gamma_2$')
                elif free=='g3' or free=='gamma_3':
                    ax.set_xlabel(r'$\gamma_3$')
                elif free=='rho' or free=='r':
                    ax.set_xlabel(r'$\rho$')
                elif free=='beta2' or free=='beta_2':
                    ax.set_xlabel(r'$\beta_2$')
                ax.set_ylabel('Difference')
    else:
        cmap = plt.cm.get_cmap('viridis', len([f'abs_diff_{p}', f'boundopt_{p}']))
        fig, axes = plt.subplots(4, 3, figsize=(16, 12), subplot_kw={'projection': '3d'})
        for row, (z, m) in enumerate(z_m_order):
            for col, p in enumerate(['P2','P3','P2P3']):
                ax = axes[row, col]
                df_filtered = df_all[(df_all['Z'] == z) & (df_all['M'] == m)]
                df_filtered = df_filtered.dropna(subset=[f'{free[0]}', f'{free[1]}', f'abs_diff_{p}', f'boundopt_{p}'])
                
                if df_filtered.empty:
                    ax.text2D(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center')
                    continue
                
                for k, nplane in enumerate([f'abs_diff_{p}', f'boundopt_{p}']):
                    plane = df_filtered
                    if plane.empty:
                        continue
                    pivot = plane.pivot(index=f'{free[1]}', columns=f'{free[0]}', values=nplane)
                    if pivot.empty:
                        continue
                    X, Y = np.meshgrid(pivot.columns.values, pivot.index.values)
                    Z = pivot.values
                    ax.plot_surface(X, Y, Z, color=cmap(k), alpha=0.45, edgecolor='none', rstride=1, cstride=1)
                
                xticks = ax.get_xticks()
                xticks_new = [tick for tick in xticks if tick in np.unique(df_filtered[f'{free[0]}'])]
                ax.set_xticks(xticks_new)
                yticks = ax.get_yticks()
                yticks_new = [tick for tick in yticks if tick in np.unique(df_filtered[f'{free[1]}'])]
                ax.set_yticks(yticks_new)
                zticks = ax.get_zticks()
                zticks_new = [round(tick,2) for tick in zticks]
                ax.set_zticks(zticks_new)
                ax.set_xlabel(f'{free[0]}', fontsize=12)
                ax.set_ylabel(f'{free[1]}', fontsize=12)
                ax.set_zlabel(f'Difference', fontsize=12)
                ax.view_init(elev=25, azim=-25)
    for i, p_name in enumerate([r'$|\psi_{P_2}^{ne} - \psi_{P_2}^{rt}$|',r'$| \psi_{P_3}^{ne} - \psi_{P_3}^{rt} |$',r'$| \psi_{P_2}^{ne} - \psi_{P_2}^{rt} | + | \psi_{P_3}^{ne} - \psi_{P_3}^{rt} |$']):
        axes[0, i].set_title(p_name, fontsize=14)
    row_titles = ["Z continuous, M continuous", "Z continuous, M binary", "Z binary, M continuous", "Z binary, M binary"]
    for row, title in enumerate(row_titles):
        fig.text(0, 0.87 - row * 0.24, title, va='center', ha='right', rotation='vertical', fontsize=12, fontweight='bold')
    plt.tight_layout()
    #plt.savefig('up/bound_across_params_2d.pdf', bbox_inches='tight')
    plt.show()
#bound_across_params_2d(conditions = {"rho":-0.2, "beta2":1.0, "g6":1.5}, free = 'a1')
#bound_across_params_2d(conditions = {"rho":-0.2, "beta2":1.0, "a1":-0.5}, free = 'g6')

#bound_across_params_2d(conditions = {"rho":-0.2, "beta2":1.5}, free=['g2','a1'])
#bound_across_params_2d(conditions = {"rho":-0.2, "beta2":1.5, "g2":1.5}, free = 'a1')
#bound_across_params_2d(conditions = {"rho":-0.2, "beta2":1.5, "a1":1.5}, free = 'g2')

#bound_across_params_2d(conditions = {"rho":-0.2, "g2":1.0}, free = ['a1', 'beta2'])
#bound_across_params_2d(conditions = {"rho":-0.2, "g2":1.0, 'beta2':1.5}, free = 'a1')
#bound_across_params_2d(conditions = {"rho":-0.2, "g2":1.0, 'a1':1.5}, free = 'beta2')

"""bound_across_params_2d(conditions = {"rho":-0.2, 'alpha_1':1.5}, free = ['gamma_3', 'beta_2'])
bound_across_params_2d(conditions = {"rho":-0.2, 'beta_2':1.5}, free = ['gamma_3', 'alpha_1'])
bound_across_params_2d(conditions = {"rho":-0.2, 'gamma_3':-1.5}, free = ['beta_2', 'alpha_1'])
bound_across_params_2d(conditions = {"rho":-0.2, 'gamma_3':-1.5, 'beta_2':1.5}, free = 'alpha_1')"""


#bound_across_params_2d(conditions = {"gamma_2":1.5, "alpha_1":-1.0, "beta_2":1.5}, free = "gamma_3")
#bound_across_params_2d(conditions = {"beta_2":1.5, "gamma_2":1.5, "gamma_3":1.5}, free = "alpha_1")

def dis_close_bound_across_param():
    df_all = pd.read_csv(f'res_test.csv')
    df_all['Z'] = df_all['mode.z']
    df_all['M'] = df_all['mode.m']
    df_all['Y'] = df_all['mode.y']
    df_all = df_all[df_all['Y']=='bi']
    df_all[f'abs_diff_P2'] = np.abs(df_all[f'diff_P2_1.0'])
    df_all[f'abs_diff_P3'] = np.abs(df_all['Ys2'] - df_all['Ys3'] - df_all['Ys2.dprime'] + df_all['Ys3.dprime'])
    for p in ['P2','P3']:        
        df_all[f'boundopt_{p}'] = df_all[[f'uppb_{p}_z_1.0',f'uppb_{p}_1_1.0',f'uppb_{p}_bo_1.0',f'uppb_{p}_bn_1.0']].min(axis=1)   
    df_all['abs_diff_P2P3'] = df_all['abs_diff_P2'] + df_all['abs_diff_P3']
    df_all['boundopt_P2P3'] = np.abs(df_all['Ys1'] - df_all['Ys3'] - df_all['Ys1.prime'] + df_all['Ys3.dprime'])
    df_all['upboundopt_P2P3'] = df_all['boundopt_P2'] + df_all['boundopt_P3']

    z_m_order = [('con', 'con'), ('con', 'bi'), ('bi', 'con'), ('bi', 'bi')]

    for sce, (z, m) in enumerate(z_m_order):
        fig, axes = plt.subplots(4, 3, figsize=(16, 12), sharey=False)
        df_filtered = df_all[(df_all['Z'] == z) & (df_all['M'] == m)]
        thresholds = []
        for col, p in enumerate(['P2','P3','P2P3']):
            df_filtered['check'] = 0
            if (p!='P2P3'):
                threshold = np.quantile(df_filtered[f'boundopt_{p}'] - df_filtered[f'abs_diff_{p}'],0.1)
                df_filtered.loc[df_filtered[f'boundopt_{p}'] - df_filtered[f'abs_diff_{p}'] <= threshold,'check'] = 1
            else:
                threshold = np.quantile(df_filtered[f'abs_diff_{p}'] - df_filtered[f'boundopt_{p}'],0.1)
                df_filtered.loc[df_filtered[f'abs_diff_{p}'] - df_filtered[f'boundopt_{p}'] <= threshold,'check'] = 1
            thresholds.append(round(threshold,3))
            for row, param in enumerate(['alpha_1','beta_2','gamma_2','gamma_3']):
                ax = axes[row,col]
                prop = (pd.crosstab(df_filtered[param], df_filtered['check'], normalize='index')
                        .reindex(columns=[0, 1], fill_value=0))
                prop.plot(
                    kind='bar',
                    stacked=True,
                    ax=ax,
                    legend=False)
                ax.set_xlabel("para")
                ax.set_ylabel("percent")
                handles, labels = ax.get_legend_handles_labels()
        row_titles = [r'$\alpha_1$',r'$\beta_2$',r'$\gamma_2$',r'$\gamma_3$']
        for row, title in enumerate(row_titles):
            fig.text(0.08, 0.8 - row * 0.2,  title, va='center', ha='right', fontsize=12, fontweight='bold')
        for i, p_name in enumerate([r'$|\psi_{P_2}^{ne} - \psi_{P_2}^{rt}$|',r'$| \psi_{P_3}^{ne} - \psi_{P_3}^{rt} |$',r'$| \psi_{P_2}^{ne} - \psi_{P_2}^{rt} | + | \psi_{P_3}^{ne} - \psi_{P_3}^{rt} |$']):
            axes[0, i].set_title(p_name + f' threshold{thresholds[i]}', fontsize=14)
        fig.legend(
            handles,
            labels,
            title="Value",
            loc="center right"
        )
        titles = ["Z continuous, M continuous", "Z continuous, M binary", "Z binary, M continuous", "Z binary, M binary"]
        fig.suptitle(titles[sce], fontsize=16)
        plt.show()

dis_close_bound_across_param()