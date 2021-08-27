# DRMS交接

## 安裝MATLAB
* 登入學校授權軟體，下載MATLAB
* https://ca.ncu.edu.tw/software/index.php
* 學校網頁有詳細安裝步驟，請按照學校提供的方式進行

## Simulation
### 實驗目錄
![](https://i.imgur.com/Y7ic01h.png)
### 操作步驟
1. main: 主程式，載入beam參數，可依照SNR、UE數量、CFO調整
    * tcp_table.mat: 不同beam適合的K
    * cfo_table.mat: 不同beam的CFO範圍
    * SNR、UE數量用array調整即可
    * CFO調整方法: 
        * main.m先改要的epsilon，圖中為2.5，輸出的檔名會帶這個參數
        * ![](https://i.imgur.com/SBra6cQ.png)
        * 改preambleFromUE.m中的function getUE()，設定UE的epsilon
        * ![](https://i.imgur.com/t6h4Pzp.png)

2. getSeqParameters: 設定sequence參數
3. getSimParameters: 設定simulation參數
4. collectPreambles: 製作不同設計的preamble pool
5. preamblesFromUE: UE選preamble、加都卜勒效應、時間延遲
    * samplePreamble: 幫UE隨機選preamble index
    * getUE: 隨機選UE要送的preamble的都卜勒效應、時間延遲
        * dopplerShift: 幫preamble加頻率偏移
    * collectSignals: preamble sequence -> preamble signal
    * passThruCh: preamble 經過通道被收端收下的訊號
    * showSelectedPreambles: 印出UE選到的preamble編號、CFO、delay
6. detectRxSignal: 偵測收到的訊號、計算時間延遲
    * detectJSAC: 偵測reference [5]
    * detectZTE: 偵測reference [8]
    * detectHSCC: 偵測DRMS
    * detect目錄: 偵測方法共用的function
        * convertL: 用root編號、cyclic shift轉成preamble編號
        * getFreqOffsetDist: 用頻偏量計算peak偏移位置
        * getPDP: Power delay profile計算
        * plotPDP: 印出來看PDP狀況，要看時第一行return註解掉
        * reversedProcess: 傅立葉轉換signal轉sequence
    * setThreshold: 設定detection threshold
7. examResult: 統計false alarm、miss detection、first access rate

### 畫圖相關
* ploting 目錄
    * auto/cross: 畫ZC的autocorrelation/cross-correlation
        ![](https://i.imgur.com/WbKaTCe.png)
    * diffCFO: 畫ZC在不同CFO下correlation變化
        ![](https://i.imgur.com/Axo7ptl.png)
    * delta_tau: 畫論文裡不同$\epsilon/K$範圍的peak偏移
        * PPP.mat: 存了delta_tau要畫的資料
        ![](https://i.imgur.com/cfIC7yk.png)
    * figureEK/figureOut/figureSNR: 畫simulation的實驗結果
        * 實驗數值請自行替換，數值取得使用parse.m
        * parse.m: 整理實驗結果，換成first access rate數值
        ![](https://i.imgur.com/asi25ya.png)
