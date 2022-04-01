/*QEM simplification(Multiple Choice Method)
  実行までの例（Eigenのパスは自身のパスを指定してください）
  g++ QEM_mcm.cpp -o QEM_mcm -std=c++17 -I /usr/local/.../include/eigen3
  ./QEM_mcm sample.off
*/
# include<iostream>
# include<fstream>
# include<cstdio>
# include<vector>
# include<map>
# include<queue>
# include<random>
# include"Eigen/Core"
# include"Eigen/Dense"
#define N 100
using Vec3d = Eigen::Vector3d;
using MatXd = Eigen::MatrixXd;
using Mat3d = Eigen::Matrix3d;

struct Face{
  int a,b,c;
  Face(int a_=0 ,int b_=0, int c_=0): a(a_), b(b_), c(c_){}
};
struct Edge{
  int v0, v1;
  Edge(int v0_=0, int v1_=0): v0(v0_), v1(v1_){}
};

struct Qv{
  Mat3d A;
  Vec3d b;
  double c;

  Qv(){
    A = Mat3d::Zero();
    b = Vec3d::Zero();
    c = 0;
  }

  Qv(const Mat3d &A_, const Vec3d & b_, const double &c_) :
    A(A_), b(b_), c(c_)
  {}

  Qv(const Vec3d &p0, const Vec3d &p1, const Vec3d &p2){
    A = Mat3d::Zero();
    b = Vec3d::Zero();
    c = 0;
    Vec3d n = (p1 - p0).cross(p2 - p0);
    double area = n.norm();
    n /= area;
    double d = -1*n.dot(p0);
    add(n, d, area);
  }

  Qv(const Vec3d &n, const double &d, const double &area){
    A << n.x()*n.x(), n.x()*n.y(), n.x()*n.z(),
         n.y()*n.x(), n.y()*n.y(), n.y()*n.z(),
         n.z()*n.x(), n.z()*n.y(), n.z()*n.z();
    A *= (area/3);
    b = (d * n) * (area/3);
    c = d * d * (area/3);
  }

  Qv operator + (const Qv &q){
    return Qv{A+q.A, b+q.b, c+q.c};
  }

  void add(const Vec3d &n, const double &d, const double &area){
    Mat3d addA;
    addA << n.x()*n.x(), n.x()*n.y(), n.x()*n.z(),
            n.y()*n.x(), n.y()*n.y(), n.y()*n.z(),
            n.z()*n.x(), n.z()*n.y(), n.z()*n.z();
    addA *= (area/3);
    A += addA;
    b += (d * n) * (area/3);
    c += d * d * (area/3);
  }

  void add(const Qv &Qf){
    A += Qf.A;
    b += Qf.b;
    c += Qf.c;
  }

  double calcQEM(const Vec3d &p){
    return p.transpose()*A*p + 2*b.dot(p) + c;
  }

  double calcQEM(){
    return -1*b.transpose()*A.inverse()*b + c;
  }
};

struct unionFind{
  //根には-1*深さ、子には親の添字を保存
  std::vector<int> node;

  unionFind(int n){
    node.resize(n, -1);
  }
  //根の添字を返す（子ノードを根につなぎ直す）
  int find(int n){
    if( node[n] < 0 ) return n;
    else return node[n] = find(node[n]);
  }

  //デバック用
  int find(int n, int cnt){
    std::cout << n << "\n";
    cnt++;
    if(cnt > 10) throw cnt;
    if( node[n] < 0 ) return n;
    else{
      node[n] = find(node[n], cnt);
      return node[n];
    }
  }

  //二つの木を深い方に統合する（深さが同じ場合xを根とする）
  void marge(int x, int y){
    x = find(x);
    y = find(y);
    int depthX, depthY;
    depthX = -1*node[x];
    depthY = -1*node[y];
    //if(depthX < 0)std::cout << depthX <<"\n";
    //if(depthY < 0)std::cout << depthY <<"\n";
    if( depthX > depthY ) node[y] = x;
    else if( depthY > depthX ) node[x] = y;
    else if( x != y ){
      node[y] = x;
      node[x]--;
    }
  }
};

int main(int argc, char *argv[]){
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

  std::vector<Vec3d> vertices;
  for(int i=0; i<v; i++){
    double x, y, z;
    char c[N];
    fgets(c, N, fl);
    sscanf(c, "%lf %lf %lf", &x, &y, &z);
    vertices.emplace_back(x,y,z);
  }

  std::vector<Face> faces;
  std::vector<Edge> edges;
  std::vector<Qv> Q;
  Q.resize(v);
  for(int i=0; i<f; i++){
    int a, x, y, z;
    char c[N];
    fgets(c, N, fl);
    sscanf(c, "%d %d %d %d", &a, &x, &y, &z);
    if( x<y ) edges.emplace_back(x, y);
    if( y<z ) edges.emplace_back(y, z);
    if( z<x ) edges.emplace_back(z, x);
    faces.emplace_back(x,y,z);

    Vec3d p0, p1, p2;
    p0 << vertices[x].x(), vertices[x].y(), vertices[x].z();
    p1 << vertices[y].x(), vertices[y].y(), vertices[y].z();
    p2 << vertices[z].x(), vertices[z].y(), vertices[z].z();

    Qv areaQf(p0, p1, p2);
    Q[x].add(areaQf);
    Q[y].add(areaQf);
    Q[z].add(areaQf);
  }

  std::map<int, double> qems;
  for(int i=0; i<edges.size(); i++){
    Qv Qtmp = Q[edges[i].v0] + Q[edges[i].v1];
    double val = Qtmp.calcQEM();
    qems.insert(std::make_pair(i,val));
  }

  std::vector<bool> stamp(v, false);
  std::vector<bool> stampe(edges.size(), false);
  unionFind u(v);
  int n = 100;
  //std::cout << "目標頂点数: "; std::cin >> n;
//QEM最小の辺と新しい頂点の算出（辺縮約をしていく）//////////////////////
  //std::random_device seed_gen;
  int seed = 123456;
  //int seed = seed_gen();
  std::mt19937 engine(seed);
  std::cout <<"seed: "<< seed <<"\n";

//ランダムに辺を複数個選ぶ
//その中でQEM最小の辺を縮約する
//辺の更新をどうするか：unionfind
//QEMの更新をどうするか：頂点に対してタイムスタンプを利用

  int dd=0;
  for(int cnt=0; cnt<v-n; cnt++){
    //std::cout <<cnt<< " roop\n";

    double min = 99999;
    int idxe = 0;
    Edge e;
    int i=0;
    while(i<5){
      int idx = engine() % edges.size();
      if( !stampe[idx] ){
      Edge tmpe = edges[idx];
      int idx0 = u.find(tmpe.v0);
      int idx1 = u.find(tmpe.v1);
      //選んだ辺のQEMが更新されてなかった場合
      if( stamp[idx0] || stamp[idx1]){
        Qv Qtmp = Q[idx0] + Q[idx1];
        double val = Qtmp.calcQEM();
        qems.at(idx) = val;
      }
      //QEMが最小か確認
      if(min > qems.at(idx)){
        min = qems.at(idx);
        e = tmpe;
        idxe = idx;
      }
      i++;
      }
      //std::cout << min <<"  ";
    }
    //std::cout <<std::endl;

    stampe[idxe] = true;
    int idx0 = u.find(e.v0);
    int idx1 = u.find(e.v1);
    stamp[e.v0] = true;
    stamp[e.v1] = true;
    Qv Qe = Q[idx0] + Q[idx1];
    Vec3d p = Qe.A.fullPivLu().solve(-1*Qe.b);
    double val = Qe.calcQEM(p);

    //辺縮約して頂点とQEMを更新
    u.marge(idx0, idx1);
    vertices[u.find(idx0)] = p;
    Q[u.find(idx0)] = Qe;
    //std::cout <<"cnt: "<< cnt << "\n";
    dd++;
  }
  std::cout << "output\n";
//メッシュの構成////////////////////////////////////////////////////////
  std::map<int,int> m;
  std::vector<Vec3d> newVertex;
  std::vector<Face> newFace;
  for(int i=0, idx=0; i<v; i++){
    if(i == u.find(i)){
      newVertex.push_back(vertices[i]);
      m.insert(std::make_pair(i, idx));
      idx++;
      dd++;
    }
  }
  std::cout <<"aa"<< dd << "\n";
  for(auto f : faces){
    int a, b, c;
    a = u.find(f.a);
    b = u.find(f.b);
    c = u.find(f.c);
    //if( a != u.find(f.a) || b != u.find(f.b) || c != u.find(f.c) ) continue;
    if( a != b && b != c && c != a){
      a = m.at(a);
      b = m.at(b);
      c = m.at(c);
      newFace.emplace_back(a, b, c);
    }
  };
//メッシュの出力//////////////////////////////////////////////////////////
  std::ofstream outFile("QEM_mcm.off");
  outFile << "OFF\n";
  outFile << newVertex.size() <<" "<< newFace.size() <<" 0\n";
  for(auto p : newVertex){
    outFile << p.x() <<" "<< p.y() <<" "<< p.z() <<"\n";
  };
  for(auto f : newFace){
    outFile << "3 " << f.a <<" "<< f.b <<" "<< f.c <<"\n";
  };

  return 0;
}
