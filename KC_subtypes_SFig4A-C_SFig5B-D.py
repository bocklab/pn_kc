
names = ["575kcy", "478kcab", "246kcprime"]

# SFig4 A-C, SFig5 B-D

for random_claws in [False, True]:
    for igroup,ana in enumerate([ana_ab, ana_prime, ana_y]):

        conn_data = ana.conn_data['glom_kc_in_claw_unit']
        ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

        if random_claws:
            stat = [get_raw_inputs(i) for i in shuffle_glom_kc_iterate(ob_conn, num_exp)]
        else:
            # random bouton null model
            stat = [get_raw_inputs(shuffle_glom_kc_w_prob(ob_conn, glom_prob)) for i in range(num_exp)]

        stat = np.array(stat)
        sd = np.nanstd(stat, axis=0)
        avg = np.nanmean(stat, axis=0)

        ob_ci = get_raw_inputs(ob_conn)
        comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)


        cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)

        gloms = df_lookup('glom_id',cm_zs.col_ids,'short_glom_name',glom_btn_table)

        reorder_idx = reorder(ClusterOrder0707, glom_idx_ids)
        t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

        fig, ax1 = plt.subplots()
        t1 = t1_zs;
        gloms = df_lookup('glom_id',t1.col_ids,'short_glom_name',glom_btn_table)
        sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, vmin=-8.53, vmax=8.53, cmap="RdBu_r")

        ax1.tick_params(bottom=False,labeltop=True, top=True, labelbottom=False)
        ax1.tick_params(axis='x',labelrotation=90)

        # run and get tbl in 191029-bouton-KC-representations_per_PN.py
        col_list = t1.col_ids
        col_colors = df_lookup('short_glom_name', gloms, 'color', tbl)

        for x in [ax1.get_xticklabels(), ax1.get_yticklabels()]:
            for idx, tick in enumerate(x):
                tick.set_color(col_colors[idx])
                if col_list[idx] in comm_ids:
                    tick.set_weight("extra bold")
        #            tick.set_bbox(dict(ec='green', fc=None, alpha=0.05))

        ax1.set_aspect("equal")
        fig.set_size_inches(16,12)
        # ax1.set_title(t8, pad=55)
        plt.show()
        # fig.savefig(save_path + t8, bbox_inches='tight')
