/*mesh simplification(vertex clustering)
  実行までの例（Eigenのパスは自身のパスを指定してください）
  g++ simplification.cpp -o simplification -std=c++17 -I /usr/local/.../include/eigen3
  ./simplification sample.off
*/
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <string>
#define N 50
/*
 * メッシュ読み込む（コマンドライン引数を利用）
 * x,y,z座標の(使わなかった：最大値)と最小値を求める（直方体でメッシュを覆う）
 * 直方体を１辺が閾値tの立方体で分ける
 * 各立方体の内部にある頂点の重心を新しい頂点にする
 * ３頂点全てが異なる立方体に含まれる面（残す面）を求める
 * メッシュの出力（simplification.off）
 */
struct Point{
  double x, y, z;
  Point(double x_=0, double y_=0, double z_=0): x(x_), y(y_), z(z_){}
  void operator += (const Point &p){
    x += p.x; y += p.y; z += p.z;
  }
  void operator /= (const double n){
    x /= n; y /= n; z /= n;
  }
};
struct Face{
  int a,b,c;
  Face(int a_=0 ,int b_=0, int c_=0): a(a_), b(b_), c(c_){}
};

//１つでも一致しているときture，全部違うときfalse
bool match(int idxX[3], int idxY[3], int idxZ[3]){
  if( idxX[0] == idxX[1] && idxY[0] == idxY[1] && idxZ[0] == idxZ[1] ) return true;
  if( idxX[1] == idxX[2] && idxY[1] == idxY[2] && idxZ[1] == idxZ[2] ) return true;
  if( idxX[2] == idxX[0] && idxY[2] == idxY[0] && idxZ[2] == idxZ[0] ) return true;
  return false;
}
int main(int argc, char *argv[]){
  //std::cout << "a\n";
//メッシュの読み込み//////////////////////////////////
  //FILE* fl = fopen("294_kitten_uniform.off", "r");
  FILE* fl = fopen(argv[1], "r");
  if(!fl){
    perror("File opennig error");
    return 0;
  }

  int v, f, tt;
  for(int i=0; i<2; i++){
    char c[N];
    fgets(c, N, fl);
    if(i == 1) sscanf(c, "%d %d %d", &v, &f, &tt);
  }

  std::vector<Point> vertex;
  for(int i=0; i<v; i++){
    double x, y, z;
    char c[N];
    fgets(c, N, fl);
    sscanf(c, "%lf %lf %lf", &x, &y, &z);
    vertex.emplace_back(x,y,z);
//      printf("%d %d %d %d\n", edge, x, y, z);
  }

  std::vector<Face> face;
  for(int i=0; i<f; i++){
    int a, x, y, z;
    char c[N];
    fgets(c, N, fl);
    sscanf(c, "%d %d %d %d", &a, &x, &y, &z);
    face.emplace_back(x,y,z);
//      printf("%d %d %d %d\n", edge, x, y, z);
  }

//（使わなかった：最大）と最小値から直方体を作る///////
  double minX(0), minY(0), minZ(0);
  for(Point &p : vertex){
    if(minX > p.x) minX = p.x;
    if(minY > p.y) minY = p.y;
    if(minZ > p.z) minZ = p.z;
  };
//各立方体の内部にある頂点の重心を新しい頂点にする////
  Point box[N][N][N];
  int b[N][N][N] = {0};
  std::vector<Point> newVertex;
  double t = 0.1;
  for(Point &p : vertex){
    int idxX, idxY, idxZ;
    idxX = (p.x - minX) / t;
    idxY = (p.y - minY) / t;
    idxZ = (p.z - minZ) / t;
    box[idxX][idxY][idxZ] += p;
    b[idxX][idxY][idxZ] += 1;
  };
  int s = 0;
  int idxBox[N][N][N] = {-1};
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      for(int k=0; k<N; k++){
        if( b[i][j][k] != 0 ){
          box[i][j][k] /= b[i][j][k];
          newVertex.push_back(box[i][j][k]);
          idxBox[i][j][k] = s;
          s++;
        }
      }
    }
  }
//３頂点全てが異なる立方体に含まれる面（残す面）を求める（面の接続関係を使うため）
  std::vector<Face> newFace;
  for(Face f: face){
    int idxX[3], idxY[3], idxZ[3];
    int fa = f.a;
    idxX[0] = (vertex[fa].x - minX) / t;
    idxY[0] = (vertex[fa].y - minY) / t;
    idxZ[0] = (vertex[fa].z - minZ) / t;
    int fb = f.b;
    idxX[1] = (vertex[fb].x - minX) / t;
    idxY[1] = (vertex[fb].y - minY) / t;
    idxZ[1] = (vertex[fb].z - minZ) / t;
    int fc = f.c;
    idxX[2] = (vertex[fc].x - minX) / t;
    idxY[2] = (vertex[fc].y - minY) / t;
    idxZ[2] = (vertex[fc].z - minZ) / t;
    if( !match(idxX, idxY, idxZ) ){
      int a = idxBox[idxX[0]][idxY[0]][idxZ[0]];
      int b = idxBox[idxX[1]][idxY[1]][idxZ[1]];
      int c = idxBox[idxX[2]][idxY[2]][idxZ[2]];
      Face ff(a,b,c);
      newFace.push_back(ff);
    }
  };
//メッシュの出力
  std::ofstream outFile("simplification.off");
  outFile << "OFF\n";
  outFile << newVertex.size() <<" "<< newFace.size() <<" 0\n";
  for(Point p : newVertex){
    outFile << p.x <<" "<< p.y <<" "<< p.z <<"\n";
  };
  for(Face f : newFace){
    outFile << "3 " << f.a <<" "<< f.b <<" "<< f.c <<"\n";
  };
  return 0;
}
