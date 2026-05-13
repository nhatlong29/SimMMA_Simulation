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
    #plt.savefig(f'rank_{sz}_Z{z}_M{m}_Y{y}{typ}.png', bbox_inches='tight')
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
    #plt.savefig(f'up/{sign_column}_{sz}_Z{z}_M{m}_Y{y}{typ}.png', bbox_inches='tight')
    plt.show()

def plot_combined(z, m, y, typ, sz, ver):
    data = pd.read_csv(f'res_Z{z}_M{m}_Y{y}{typ}.csv')
    data['tol_re'] = data[[f'P1_{ver}',f'P2_{ver}',f'P3_{ver}',f'P4_{ver}',f'Ze_{ver}']].sum(axis=1)
    data['per_Ze'] = np.abs(data[f'Ze_{ver}'])/np.abs(data['tol_re'])
    data['mag_re'] = np.sum(np.abs(data[[f'P1_{ver}',f'P2_{ver}',f'P3_{ver}',f'P4_{ver}',f'Ze_{ver}']]),axis=1)
    data['per_magZe'] = np.abs(data[f'Ze_{ver}'])/np.abs(data['mag_re'])
    if typ=='':
        r = pd.read_csv(f'r_Z{z}_M{m}_Y{y}{typ}.csv')
        r = pd.merge(r,data[['per_Ze','per_magZe']],left_index=True,right_index=True)
        r[f'Ze_{ver}'] = np.abs(r[f'Ze_{ver}'])
    else:
        r = pd.read_csv(f'r_Z{z}_M{m}_Y{y}{typ}.csv')
        r = pd.merge(r,data[['per_Ze','per_magZe',f'diff_P2_{ver}']],left_index=True,right_index=True)
        r[f'Ze_{ver}'] = np.abs(r[f'Ze_{ver}'])
        r[f'diff_P2_{ver}'] = np.abs(r[f'diff_P2_{ver}'])
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
    data[f'diff_P2_{ver}'] = np.abs(data[f'diff_P2_{ver}'])
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
    plt.savefig(f'up/{sz}_Z{z}_M{m}_Y{y}{typ}_ver{ver}.png', bbox_inches='tight')
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
"""
sz = None
df_all = pd.DataFrame()
for z in ['bi','con']:
    for m in ['bi','con']:
        for y in ['bi','con']:
            for ver in ['1.0','2.0','3.0']:
                df = plot_combined(z=z, m=m, y=y, sz=sz, ver=ver)
                #df_all = pd.concat([df_all,df], axis=0, ignore_index=True)
"""
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
    df_all['psiz'] = np.abs(df_all['bound_P2_1.0'] - np.abs(df_all['diff_P2_1.0']))
    df_all['probb'] = np.abs(df_all['bbound_P2_1.0'] - np.abs(df_all['diff_P2_1.0']))
    df_all.loc[df_all['psiz']>=df_all['probb'],'check'] = 1
    df_all.loc[df_all['psiz']<df_all['probb'],'check'] = 0
    
    for z in ['bi','con']:
        for m in ['bi','con']:
            print(f'z {z} m {m}: {df_all.loc[((df_all['Z']==z)&(df_all['M']==m)),'check'].mean()}')
    
    cols = ['abs_diff', 'psiz', 'probb']
    df_all["group"] = ("Z" + df_all["Z"].astype(str) +"_M" + df_all["M"].astype(str))
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
    fig, axes = plt.subplots(2, len(groups), figsize=(18, 8), sharey=False)
    thr = df_all["abs_diff"].median()
    # -------------------------
    # ROW 1: abs_diff → probb
    # -------------------------
    for j, grp in enumerate(groups):
        ax = axes[0, j]
        sub = df_long[df_long["group"] == grp]
        wide = sub.pivot(index="id", columns="Variable", values="Value")
        abs_vals = wide["abs_diff"]
        for i, row in wide.iterrows():
            color = "red" if abs_vals.loc[i] < thr else "blue"
            ax.plot(
                ["abs_diff", "probb"],
                [row["abs_diff"], row["probb"]],
                color=color,
                alpha=0.2,
                linewidth=1
            )
        ax.scatter(wide.index.map(lambda _: "abs_diff"),
                wide["abs_diff"], s=10)
        ax.scatter(wide.index.map(lambda _: "probb"),
                wide["probb"], s=10)
        ax.set_title(grp)

    # -------------------------
    # ROW 2: abs_diff → psiz
    # -------------------------
    for j, grp in enumerate(groups):
        ax = axes[1, j]
        sub = df_long[df_long["group"] == grp]
        wide = sub.pivot(index="id", columns="Variable", values="Value")
        abs_vals = wide["abs_diff"]
        for i, row in wide.iterrows():
            color = "red" if abs_vals.loc[i] < thr else "blue"
            ax.plot(
                ["abs_diff", "psiz"],
                [row["abs_diff"], row["psiz"]],
                color=color,
                alpha=0.2,
                linewidth=1
            )
        ax.scatter(wide.index.map(lambda _: "abs_diff"),
                wide["abs_diff"], s=10)
        ax.scatter(wide.index.map(lambda _: "psiz"),
                wide["psiz"], s=10)
    axes[0,0].set_ylabel("abs_diff → probb")
    axes[1,0].set_ylabel("abs_diff → psiz")
    plt.tight_layout()
    plt.show()

    """
    plt.figure(figsize=(10,6))
    for _, subdf in df_long.groupby("id"):
        plt.plot(
            subdf["Variable"],
            subdf["Value"],
            alpha=0.2,
            linewidth=1
        )
    sns.scatterplot(
        data=df_long,
        x="Variable",
        y="Value",
        hue="group",
        s=40
    )
    plt.title("Connected variables within rows")
    plt.tight_layout()
    plt.show()"""

count_bound()

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
    #plt.savefig(f'up/a0_a1_ ranging.png', bbox_inches='tight')
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
    plt.savefig(f'up/bound{typ}.png', bbox_inches='tight')
    #plt.show()

#diff_p2_range_ver1(typ=2)