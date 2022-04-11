## Analyzing gain/loss of AMR genes (title to be determined)

#### TOC

- [NK_A0021](#nk_a0021) : *AMR* Future prediction of ARG dissemination
- [NK_A0022](#nk_a0022) : *AMR* KOfamScan for downloaded genomes
- [NK_A0023](#nk_a0023) : *Cross domain* Feature selection by RF after rough selection by MI
- [NK_A0024](#nk_a0024) : *Cross domain* Random ForestのTreeのdepthとAUCの関係性を見てみる
- [NK_A0025](#nk_a0025) : *AMR* Random ForestのTreeのdepthとAUCの関係性を見てみる AMR
- [NK_A0026](#nk_a0026) : *AMR* 特徴量選択を行なって特徴量個数とAUCの関係性を見てみる AMR
- [NK_A0027](#nk_a0027) : *AMR* 昔の枝と最近の枝で進化パターンが異なるのかどうかを検証する
- [NK_A0028](#nk_a0028) : *AdhE* apillactobacilus, fructobacillusの系統樹抽出
- [NK_A0029](#nk_a0029) : *AMR* ARGのdN/dS解析
- [NK_A0030](#nk_a0030) : *Cross domain* Visualize phylogenetic distribution of genes
- [NK_A0031](#nk_a0031) : *Cross domain* Synteny analysis **Including raw data to be backupped**
- [NK_A0032](#nk_a0032) : *AMR* Same experiment as A0027 with randomized branch division
- [NK_A0033](#nk_a0033) : *Cross-domain* Aerobic/anaerobicな環境で有利と思われる遺伝子を同定する
- [NK_A0034](#nk_a0034) : *AMR* 未来予測やり直し（特徴量選択の枠組みなど最適化したので）
- [NK_A0035](#nk_a0035) : *AMR* 未来予測やり直し（過去の進化と最近の進化に分けて学習）
- [NK_A0036](#nk_a0036) : *AdhE* Fructobacillaceaeの系統樹上での祖先状態推定
- [NK_A0037](#nk_a0037) : *Cross domain* AtoA, BtoB, AtoB, BtoAの特徴量選択なしの進化予測。NK_A0015のやり直し
- [NK_A0038](#nk_a0038) : *AdhE* Diversitreeの練習
- [NK_A0039](#nk_a0039) : *AdhE* Diversitreeの本番
- [NK_A0040](#nk_a0040) : *nosZ* 按田さん解析 nosZ系統分布の可視化
- [Memo](#memo)

### NK_A0021

- **Objective**

    - 各種が各ARGを将来獲得する確率を予測する

- **Command**

    - **Computer** Crux

    - (*2021/10/26*) 未来予測

        - https://github.com/IwasakiLab/PredictingAMRevolution / [55c4a91](https://github.com/IwasakiLab/PredictingAMRevolution/commit/55c4a91d60bebca9f3a2291d0fecbf32f1d9d4dd)

        - https://github.com/UTNK/handycsv / [c455d09](https://github.com/UTNK/handycsv/commit/c455d09be7f9b9cbccbd9f1d72b15031e336bce8)


        ```shell
        (base) apfe1{naoki}1028: qsub -J 0-249 '/user1/scl9/naoki/aptmp_naoki/code/PredictingAMRevolution/shell/NK_A0021_future_prediction.amr.sh'
        1243328[].apfe3
        ```

    - (*2021/10/27*) 集計

        ```shell
        (base) apfe1{naoki}1003: cat future/abko/bacteria/mlgtdb_MPPA/*/gain/sp_prob.*_gain_RF_RF.*.test.txt > results/sp_prob.gain_RF_RF.txt &
        [1] 23320
        (base) apfe1{naoki}1010: gzip results/sp_prob.gain_RF_RF.txt
        ```


### NK_A0022

- **Objective**

    - ダウンロードした同種別ゲノムに対してKOfamScanを行う

- **Command**

    - (*2021/10/21*) KOfamScanのインストール

        ```shell
        (amr) apfe1{naoki}1058: wget https://www.genome.jp/ftp/tools/kofam_scan/kofam_scan-1.3.0.tar.gz
        (base) apfe1{naoki}1023: cd kofam_scan/
        (base) apfe1{naoki}1069: mkdir download
        (base) apfe1{naoki}1070: cd download/
        (base) apfe1{naoki}1071: wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
        (base) apfe1{naoki}1072: wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
        (base) apfe1{naoki}1076: gunzip ko_list.gz
        (base) apfe1{naoki}1077: tar xvf profiles.tar.gz
        ```

    - (*2021/10/21*) より多くのゲノムをダウンロードできたダウンロード方法2: 直接URLを指定してダウンロードによってダウンロードしたゲノムに対してKOfamScanを行った

        - https://github.com/IwasakiLab/PredictingAMRevolution / [e18c9ef](https://github.com/IwasakiLab/PredictingAMRevolution/commit/e18c9ef9c1671b59f32ac06d3b6bb7a5f0bac52a)

        ```shell
        (base) apfe1{naoki}1003: qsub -J 0-999 '/user1/scl9/naoki/aptmp_naoki/code/PredictingAMRevolution/shell/NK_A0022_run_KOfamScan.sh'
        1237974[].apfe3
        ```

    - (*2021/10/22*) KOfamScanの結果をまとめた

        ```shell
        (base) apfe1{naoki}1008: for dir in $(ls result/); do cat result/${dir}/kofamscan.* | gunzip | grep '\*' | cut -f 2,3 | awk -v dir=$dir '{print dir"\t"$0}'; done | gzip > concat/concat.result.txt.gz &
        [1] 25770
        ```

    - **Computer** Local

        ```shell
        (base) konnonaoki@KonnonoMacBook-Pro ~ % scp -r3 kyodai:/user1/scl9/naoki/aptmp_naoki/data/PredictingAMRevolution/NK_A0022/concat crux:/data/naoki-konno/data/PredictingAMRevolution/NK_A0022
        ```

    - **Computer** Crux

        - https://github.com/UTNK/handytree / [928cf86](https://github.com/UTNK/handytree/commit/928cf86d177420c7c694ebcdc3bd446e0a7e1aa7)

        ```shell
        (base) [naoki-konno@crux NK_A0022]$ less tables/NK_A0022_possession_ratio_pred.txt | cut -f2 | tail -n +2| sort | uniq > tables/analyzed_species.txt
        (base) [naoki-konno@crux NK_A0022]$ less tables/analyzed_species.txt |wc -l
        1125
        (base) [naoki-konno@crux NK_A0022]$ mkdir tree
        (base) [naoki-konno@crux NK_A0022]$ tree_extract -t ../NK_A0004/bac120_r202.selected.internal_renamed.tree -n tables/analyzed_species.txt > tree/bac120_r202.selected.internal_renamed.pangenome.tree &
        [1] 21652
        ```

    - (*2021/11/07*) Accession IDとspecies nameを紐付けた後、パンゲノム解析の対象種の中からPATRICに登録された重要な病原性種 (NK_A0020参照)を絞り込んだ

        - https://github.com/UTNK/handytree / [928cf86](https://github.com/UTNK/handytree/commit/928cf86d177420c7c694ebcdc3bd446e0a7e1aa7)

        ```shell
        (base) [naoki-konno@crux NK_A0022]$ cat ../NK_A0001/bac120_metadata_r202.tsv | cut -f1,17,79 | tr ';' '\t'|cut -f1,8,15|sed 's/s__//g' | tail -n +2> tables/acc_gtdbsp_ncbisp.txt
        (base) [naoki-konno@crux NK_A0022]$ csvjoin -c1,3 --no-header-row -t ../NK_A0020/tables/pathogen_list.txt tables/acc_gtdbsp_ncbisp.txt | tail -n +2 | awk -F, '{print $2"\t"$1"\t"$3}'> tables/acc_gtdbsp_ncbisp.pathogen.txt
        (base) [naoki-konno@crux NK_A0022]$ csvjoin -c1,1 --no-header-row -t tables/analyzed_species.txt tables/acc_gtdbsp_ncbisp.pathogen.txt | tail -n +2 | tr ',' '\t' > tables/acc_gtdbsp_ncbisp.pathogen.analyzed.txt
        (base) [naoki-konno@crux NK_A0022]$ less acc_gtdbsp_ncbisp.pathogen.analyzed.txt | wc -l
        186
        (base) [naoki-konno@crux NK_A0022]$ less acc_gtdbsp_ncbisp.pathogen.analyzed.txt | cut -f2 | sort |uniq | wc -l
        121
        ```

    - (*2021/11/07*) 絞り込んだ種の系統樹を抜き出した

        ```shell
        (base) [naoki-konno@crux NK_A0022]$ cat tables/acc_gtdbsp_ncbisp.pathogen.analyzed.txt | cut -f1 > tables/analyzed_species.pathogen.txt
        (base) [naoki-konno@crux NK_A0022]$ tree_extract -t tree/bac120_r202.selected.internal_renamed.pangenome.tree -n tables/analyzed_species.pathogen.txt > tree/bac120_r202.selected.internal_renamed.pangenome.pathogen.tree
        before extraction: 1125 tips
        after extraction: 186 tips
        ```