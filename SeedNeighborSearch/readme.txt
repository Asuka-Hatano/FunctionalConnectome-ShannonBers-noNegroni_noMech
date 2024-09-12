Copied from local folder:
Trial014_RecalcTau_SeedNeighborSearch

Main folder (upper folder) is the main stream, most recent code.

Only the SeedNeighbor Search section of the program should be adopted, and when reusing, it should be merged accordingly with the main stream.

* Tau calculation has been corrected
* The convergence condition remains outdated (it is not using the correct Fractional Error).
* Main folder output includes additional columns, number of columns in the output data is different.

The calculations in this folder output the results as Yfinal and Metrics files each time.
While the files are numbered, there are missing numbers. This is due to the use of parfor, the presence of randomness in SeedNeighbor, and occasional crashes during Tau calculation.

CorrectMetricsYfinals***.m code consolidates the results for only the existing numbered files into a single file.

--------------
Tauの計算は修正されているものの、収束条件が古いまま（正しいFractional Errorになっていない）
本フォルダの方のプログラムでは追加の列も計算していて、出力データの列数も異なる。

SeedNeighbor Searchの部分のプログラムのみを採用して、再度利用する際にはマージすべし。

本フォルダの計算は、結果のYfinal と Metricsファイルを毎回出力する。
ナンバリングしているが、抜け番号がある。（Parforで呼ぶのと、SeedNeighborには乱数があり全部の番号が実行されないのと、Tau計算のせいなのか稀に計算が落ちるので）

CorrectMetricsYfinals のコードは、通し番号のうち、存在する番号のみの結果を一つのファイルにまとめる。小さなテキストファイルが大量に存在するのは扱いづらく、読み込みも面倒なので必要。
