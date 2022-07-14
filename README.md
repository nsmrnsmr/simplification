# mesh simplification
メッシュの簡略化を行うプログラム

入力：　メッシュ（off形式），縮約する辺の本数

出力：　簡略化したメッシュ（off形式）

### 実行方法
- g++ simplification.cpp -o simplification -std=c++17 -I /usr/local/.../include/eigen3
- ./simplification sample.off
  
---
<p align="center">
  <img src="image/mesh00.png" width="">
  <br>
  <em>元となるメッシュ</em>
  <img src="image/vertexClustering00.png" width="">
  <br>
  <em>頂点クラスタリングを用いた簡略化の結果</em>
  <img src="image/QEMsimplification100.png" width="">
  <br>
  <em>QEMを用いた簡略化の結果</em>
</p>
