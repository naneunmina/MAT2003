#include <iostream>
#include <random>
#include <cmath>

using namespace std;

//행렬 클래스
class matrix{
   public:
   int row, col;
   double mat[8][8];

   matrix(int row, int col);                               //constructor
   matrix rowvec(int num);          //행렬에서 num행의 벡터 한개만 가져옴
   matrix operator*(double a);                            //행렬에 상수곱
   double operator+(const matrix inner);                           //내적
   void rannum();                           //행렬을 랜덤한 값으로 채워줌
   void printmat();                                           //행렬 출력
   double vectorsize(int a);                              //벡터크기 반환
   bool is_independent(); //gauss elimination사용. pivot 0일시 false 반환
   void GramSchmidt();                                    //Gram-Schmidt
   void plus(matrix b);                                       //행렬+행렬
   void mi(int r, matrix b);                                  //행렬-행렬
   void makevec(matrix b);                          //linear combination
   matrix trans();                                                 //전치
};

//main
int main(){
   int n;
   cout << "n : ";
   cin >> n;

   matrix stan(n,n),ut(n,n);
   while (1){
      stan.rannum();
      ut=stan;
      if (ut.is_independent()) break;
   }
   cout << "\n" << n << " Vector :\n";
   stan.printmat();
   cout << "\n" << " Upper Triangle :\n";
   ut.printmat();
   matrix gram=stan;
   gram.GramSchmidt();
   cout << "\n" << " Gram-Schmidt :\n";
   gram.printmat();
   cout << "\n";

   for (int i=1; i<=n; i++){
      cout << i << "번째 vector의 자기자신 내적: " << gram.rowvec(i)+gram.rowvec(i).trans() << "\n";
   }
   cout << "\n";
   for (int i=1; i<=n-1; i++){
      for (int j=i+1; j<=n; j++){
         cout << i << "번째 vector와 " << j << "번째 vector의내적: " << gram.rowvec(i)+gram.rowvec(j).trans() << "   반올림: " << round((gram.rowvec(i)+gram.rowvec(j).trans())*100)/100 << "\n";
      }
   }
   for (int i=0; i<3; i++){
      cout << "\n";
      matrix proof(1,n);
      proof.rannum();
      gram.makevec(proof);
   }
   return 0;
}


//행렬 클래스 정의
matrix::matrix(int row, int col){
   this->row =row;
   this->col=col;
}
matrix matrix::operator*(double a){
   for (int i=1; i<=this->col; i++){
      this->mat[1][i]*=a;
   }
   return *this;
}
double matrix::operator+(const matrix inner){
   double ret=0;
   for (int i=1; i<=inner.row; i++) ret+=(this->mat[1][i])*(inner.mat[i][1]);
   return ret;
}
matrix matrix::rowvec(int num){
   matrix ret(1,this->col);
   for (int i=1; i<=this->col; i++){
      ret.mat[1][i]=this->mat[num][i];
   }
   return ret;
}
void matrix::rannum(){
   random_device rd;
   mt19937 ran(rd());
   uniform_int_distribution<int> dis(0,9);
   for (int i=1; i<=this->row; i++){
      for (int j=1; j<=this->col; j++) this->mat[i][j]=dis(ran);
   }
}
void matrix::printmat(){
   for (int i=1; i<=this->row; i++){
      for (int j=1; j<=this->col; j++){
      cout << this->mat[i][j] << " ";
      }
      cout << "\n";
   }
}
double matrix::vectorsize(int a){
   double ret;
   for (int i=1; i<=this->col; i++) ret+=this->mat[a][i]*this->mat[a][i];
   return sqrt(ret);
}
bool matrix::is_independent(){
   for (int i=1; i<=row; i++){
      for (int j=i+1; j<=col; j++){
         double l=mat[j][i]/mat[i][i];
         for (int k=i; k<=col; k++){
            mat[j][k]-=mat[i][k]*l;
         }
      }
      if(mat[i][i]==0) return false;
   }
   return true;
}
void matrix::GramSchmidt(){
   for (int i=1; i<=row; i++){
      for (int j=1; j<i; j++){
         this->mi(i,rowvec(j)*(rowvec(i)+rowvec(j).trans()));
      }
      double si=vectorsize(i);
      for (int j=1; j<=col; j++) mat[i][j]/=si;
   }
}
void matrix::plus(matrix b){
   for (int i=1; i<=col; i++){
      mat[1][i]+=b.mat[1][i];
   }
}
void matrix::mi(int r, matrix b){
   for (int i=1; i<=col; i++){
      mat[r][i]-=b.mat[1][i];
   }
}
void matrix::makevec(matrix b){
   cout << "randomly generated vector\n";
   b.printmat();
   matrix ans(1,col);
   for (int i=1; i<=col; i++) ans.mat[1][i]=0;
      for (int i=1; i<=row; i++){
         cout << i << "번째 xe : ";
         (rowvec(i)*(b+rowvec(i).trans())).printmat();
         ans.plus(rowvec(i)*(b+rowvec(i).trans()));
      }
      cout << "모든 xe의 합 : ";
      ans.printmat();
      cout << "\n";
}
matrix matrix::trans(){
   matrix ret(col,row);
      for (int i=1; i<=row; i++){
         for (int j=1; j<=col; j++) ret.mat[j][i]=mat[i][j];
      }
   return ret;
}