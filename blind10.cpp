#ifndef TANK2_RNG_H
#define TANK2_RNG_H
#include <stdint.h>
#include <random>
#include <limits>
#include <cstring>
namespace RNG {
	using namespace std;
	typedef uint64_t u64;
	typedef uint32_t u32;
	typedef int64_t i64;
	typedef int32_t i32;
	struct xss1k {
		typedef u64 result_type;
		static const u64 rand_max = ~u64(0);
		u64 s[16], p;
		xss1k() :
			s{0x79a3d197a3eda564ull, 0xea0832f5f13eca6eull,
			  0x9b25d666c98ad37eull, 0xc47e8b052d0159dfull,
			  0xc3d03d041c267d8dull, 0xfea949190d0bf9f8ull,
			  0x07f1e5bc39068da5ull, 0x44679e2a74432b86ull,
			  0x205ca0fd8db63493ull, 0x3c086aaa6dd55a74ull,
			  0xc4e064b3e570282bull, 0x4c7b4621df4da54dull,
			  0x7b14f5648cbe10cdull, 0x848656ce83c494a5ull,
			  0xb02648187107686dull, 0x964d7d7714d55e34ull},
			p(0) {
		}// a really long initializer list; in fact generated with /dev/random
		template<typename T>
		xss1k(const T *buffer_type) { memcpy(s, buffer_type, sizeof s); }
		template<typename T>
		inline void seed(const T *buffer_type) { memcpy(s, buffer_type, sizeof s); }
		inline u64 operator()() {
			register u64 s0 = s[p], s1 = s[++p &= 15];
			s1 ^= s1 << 31;
			return (s[p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30)) * 0x106689d45497fdb5ull;
		}
	};
	struct xsadw {
		typedef u64 result_type;
		static const u64 rand_max = ~u64(0);
		u64 s[2];
		xsadw() : s{0xc4f7ba40435d8ee7ull, 0x54d21665a4290c80ull} {}
		template<typename T>
		xsadw(const T *buffer_type) { memcpy(s, buffer_type, sizeof s); }
		template<typename T>
		inline void seed(const T *buffer_type) { memcpy(s, buffer_type, sizeof s); }
		inline u64 operator()() {
			u64 x = *s, y = s[1];
			x ^= x << 23;
			return (s[1] = x ^ (*s = y) ^ (x >> 17) ^ (y >> 26)) + y;
		}
	};
	typedef unsigned char u8;
	template<class T>
	struct Pool {
		T *rng;
		int noko;
		using result_type=typename T::result_type;
		static const size_t length = 16 * sizeof(result_type);
		union {
			struct { typename T::result_type t[16]; } dat;
			struct { u8 t[length]; } buf;
		} b;
		Pool(T *a) : rng(a), noko(length) {}
		void regen() {
			noko = 0;
			for (int i = 0; i < 16; ++i)
				b.dat.t[i] = rng->operator()();
		}
		u8 operator()() {
			if (noko == length)regen();
			return b.buf.t[noko++];
		}
		void operator()(u8 *a, u8 *b) {
			while (a != b)
				*a++ = this->operator()();
		}
		u64 extrude_bit(int number) {
			if (number <= 0)return 0;
			u64 mask = number >= 64 ? u64(-1) : ((1u << u64(number)) - 1);
			u64 ans = 0;
			while (number >= 0) {
				ans = ans << 8 | operator()();
				number -= 8;
			}
			return ans & mask;
		}
	};
	RNG::xsadw rnd;
	RNG::Pool<RNG::xsadw> *rnpool;
	template<class Rng_Type>
	struct unsigned_uniform_dist {
		using result_type=typename Rng_Type::result_type;
		static const result_type tmax = numeric_limits<result_type>::max();
		result_type usable_max, upper_limit;
		Rng_Type *rng;
		unsigned_uniform_dist(Rng_Type *rng, result_type ul) : rng(rng), upper_limit(ul) {
			if (ul) {
				usable_max = (tmax / ul) * ul - 1;
				if ((usable_max + ul) > usable_max)
					usable_max += ul;
			}
		}
		result_type operator()() {
			result_type c = rng->operator()();
			if (!upper_limit)return c;
			return c <= usable_max ? c % upper_limit : operator()();
		}
	};
	template<class rng>
	auto unsigned_sample(rng *a, typename rng::result_type ul) {
		return (unsigned_uniform_dist<rng>(a, ul))();
	}
	template<class t>
	int select(int len, t a, t h) {
		int d = 0;
		for (int i = 0; i < len; ++i)
			d += a[i];
		int sample = unsigned_sample(&rnd, u64(d));
		for (int i = 0; i < len; ++i)
			if ((sample -= a[i]) < 0)
				return h[i];
	}
	void srandit() {
		u32 r[4];
		std::random_device rndv;
		r[0] = rndv();
		r[1] = rndv();
		r[2] = rndv();
		r[3] = rndv();
		rnd.seed(r);
		rnpool = new RNG::Pool<RNG::xsadw>(&rnd);
	}
}
#endif //TANK2_RNG_H
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<algorithm>
#include<vector>
#include<random>
#include"jsoncpp/json.h"
#include<chrono>
#include <iostream>
using namespace std;
using namespace RNG;
namespace Complete_Simulator {
	using std::swap;
	typedef unsigned char u8;
	const u8
		space = 0x00,
		water = 0x01,
		brick = 0x02,
		steel = 0x03,
		headq = 0x04,
		tank00 = 0x08,
		tank01 = 0x10,
		tank10 = 0x20,
		tank11 = 0x40,
		tanks[4] = {0x08, 0x10, 0x20, 0x40},
		typemask = 0x03,
		stopping = 0x7e,
		tankmask = 0x78;
	const int dx[4] = {0, 1, 0, -1}, dy[4] = {-1, 0, 1, 0};
	const int exdis[8] = {0, 255, 1, 255, 0, 255, 255, 255};
	inline bool in_range(int a) { return a >= 0 && a <= 8; }
	inline bool move_move(int a) { return a < 4; }
	inline bool mono(u8 a) {
		a = a & 0xf8;
		return (a & (a - 1)) == 0;
	}
	const u32 LTank = 0x30210;
	inline int lowbit(u8 t) { return t & (-t); }
	inline int TankId(u8 t) { return (LTank >> (t >> 2)) & 3; }
	const u8 inilis[4] = {32u, 96u, 104u, 40u};
	struct pos { u8 x:4, y:4; };
	ostream &operator<<(ostream &a, pos b) { return a << "(" << (int) b.x << "," << int(b.y) << ")"; }
	const pos dire[4] = {{0,   0xf},
								{1,   0},
								{0,   1},
								{0xf, 0}};
	const pos home[2] = {{4, 0},
								{4, 8}};
	bool valid(pos t) { return t.x < 9 && t.y < 9; }
	inline pos operator+(pos a, pos b) { return (pos) {u8(a.x + b.x), u8(a.y + b.y)}; }
	inline pos &operator+=(pos &a, pos b) {
		a.x += b.x;
		a.y += b.y;
		return a;
	}
	inline int maha(pos a, pos b) { return abs((int) a.x - (int) b.x) + abs((int) a.y - (int) b.y); }
	inline bool operator==(pos a, pos b) { return a.x == b.x && a.y == b.y; }
	const char *out[9] = {"?¡è", "?¨¹", "?¨²", "?y", "??", "?", "?", "?", "?"};
	/**
	 * Vector infrastructure notation:
	 * board::u8 as positioning : [Hi][7-4 x][3-0 y][Lo]
	 * pos : [Hi][7-4 y][3-0 x][Lo] (bitfield structure)
	 * main::Dist and for all serialized 2d matrix access: y*9+x
	 * pos are recommended usage as it is high-level API and
	 * has convenient operation semantics.
	 */
	pos makepos(int x, int y) { return (pos) {u8(x), u8(y)}; }
	inline pos movedir(pos a, int move) { return ~move ? move > 3 ? a : a + dire[move] : a; }
	struct board {
		/**
		 * board.dick: [Hi][7-4 Not dead?][3-0 Can't shoot?][Lo]
		 */
		u8 m[9][9], t[4], dick;
		u8 &operator[](pos t) { return m[t.y][t.x]; }
		u8 operator()(int x, int y) { return m[y][x]; }
		u8 &operator[](u8 g) { return m[g & 15][g >> 4]; }
		bool ise(int tnok, int tn2 = 0) {
			if (tnok > 1)return tank(tnok).y < 4 + tn2;
			else
				return tank(tnok).y > 4 - tn2;
		}
		pos shopo(pos t, int dir) {
			for (; valid(t += dire[dir]);)
				if (operator[](t) & stopping)
					return t;
			return makepos(4, 4);
		}
		int relY(int tnk) { return tnk & 2 ? t[tnk] & 15 : (8 - (t[tnk] & 15)); }
		bool canShoot(pos b) { return (~dick) & (operator[](b) >> 3); }
		bool safe(pos a, pos b) {
			if (a.x != b.x && a.y != b.y)return true;
			if (shopo(a, shotTo(a, b)) == b && canShoot(b))
				return false;
			return true;
		}
		bool csafe(pos a, pos b) {
			if (a.x != b.x && a.y != b.y)return true;
			if (shopo(a, shotTo(a, b)) == b)
				return false;
			return true;
		}
		int shotTo(pos a, pos b) {
			if (a.x < b.x)return 1;
			if (a.x > b.x)return 3;
			if (a.y < b.y)return 2;
			if (a.y > b.y)return 0;
			return -1;
		}
		int TypeOfSemi(int tankMe, int move, int tankOpp) {
			const pos thi = tank(tankMe);
			const pos tha = tank(tankOpp);
			if (thi == tha)return 3;
			return maha(thi, tha) == 2 ? 1 : 2;
		}
		int dangerous(int tankMe, int move, int tankOpp) {
			/// safe: after taking the move, the tank won't be shot (0)
			/// conditional safe: after taking the move, the position could be shot in the orthogonal direction (1)
			/// dangerous: after taking the move, the position will definitely get shot [shot in the parallel](2)
			if (!canMove(tankOpp))return 0;
			const pos me = tank(tankMe), opp = tank(tankOpp), to = movedir(me, move);
			if (me == opp) {
				if (move == -1)return 0;
				if (move > 3)return 0;
				return 1;
			}
			if (move == -1) return safe(me, opp) ? 0 : 2;
			if (move < 4) {
				/// a move step
				if (!canShot(tankOpp))return 0;
				operator[](me) ^= tanks[tankMe];
				if ((shopo(to, move ^ 1) == opp) || (shopo(to, move ^ 3) == opp))
					return (operator[](me) ^= tanks[tankMe]), 1;
				if ((shopo(to, move) == opp) || (shopo(to, move ^ 2) == opp)) return (operator[](me) ^= tanks[tankMe]), 2;
				operator[](me) ^= tanks[tankMe];
				return 0;
			} else {
				/// a shoot step
				if (canShot(tankOpp))
					for (int i = 1; i < 4; ++i)
						if (shopo(me, (move ^ i) - 4) == opp)
							return 2;
				/// If the tank could shoot me after I took the shot
				/// Only happens when I shot a brick
				const pos sh = shopo(me, move - 4);
				//	cerr << me<< sh << endl;
				if (operator[](sh) == brick) {
					/// remove the brick
					operator[](sh) = 0;
					/// in this case, opponent tank can only move in from
					/// orthogonal direction, and shoot me in my front
					const pos m1 = movedir(opp, move ^ 5);
					const pos m2 = movedir(opp, move ^ 7);
					if (shopo(opp, move ^ 6) == me || ((!operator[](m1)) && (shopo(m1, move ^ 6) == me)) ||
						 ((!operator[](m2)) && (shopo(m2, move ^ 6) == me))) {
						//		cerr << "Here" << endl;
						/// If I can't hide, then I am in danger.
						/// I can hide only when I can move to a block which
						/// Another opponent tank can't shoot
						if (((tankMe&2)&&(me.y<4)) ||(me.y>4)&&(!(tankMe&2))) {
							const pos to1 = movedir(me, move ^ 5), to2 = movedir(me, move ^ 7), oppt = tank(tankOpp ^ 1);
							if ((operator[](to1)) && (operator[](to2)))goto _l_dang;
							if (!canMove(tankOpp ^ 1))goto _l_safe;
							if (!(operator[](to1)) && csafe(to1, oppt))goto _l_safe;
							if (!(operator[](to2)) && csafe(to2, oppt))goto _l_safe;
						}
						goto _l_dang;
					}
					/// Re-add the brick
				_l_safe:
					operator[](sh) = brick;
					return 0;
				_l_dang:
					operator[](sh) = brick;
					return 2;
				}
				return 0;
			}
		}
		int dangerous(int tankMe, int move) {
			int c1 = dangerous(tankMe, move, tankMe ^ 2);
			int c2 = dangerous(tankMe, move, tankMe ^ 3);
			//	cerr << move << " is " << c1 << "," << c2 << endl;
			if (c1 == 2 || c2 == 2)return 2;
			if (c1 == 0 && c2 == 0)return 0;
			return 1;
		}
		bool contain(pos a, pos b, u8 mask) {
			const int d = shotTo(a, b);
			for (; !(a == b);)
				if (operator[](a += dire[d]) & mask)
					return true;
			return false;
		}
		bool containSteel(pos a, pos b) {
			const int d = shotTo(a, b);
			for (; !(a == b);)
				if (operator[](a += dire[d]) == steel)
					return true;
			return false;
		}
		bool canGoLine(int y, int x, int xt) { return !contain(makepos(x, y), makepos(xt, y), 1 | tankmask); }
		bool canMoveLine(int y, int x, int xt) { return !contain(makepos(x, y), makepos(xt, y), 0xff); }
		int shoodire(pos a, pos b) {
			for (int i = 0; i < 4; ++i)
				if (shopo(a, i) == b)
					return i ^ 2;
			return -1;
		}
		pos tank(int i) { return (pos) {u8(t[i] >> 4), u8(t[i] & 15)}; }
		inline bool do_move(int a, int v, u8 om[9][9]) {
			if (v == -1)return false;
			int x = t[a] >> 4, y = t[a] & 15;
			int Dx = x + dx[v], Dy = y + dy[v];
			if (in_range(Dx) && in_range(Dy) && ((typemask & m[Dy][Dx]) == space)) {
				if (om[Dy][Dx])return true;
				m[y][x] ^= tanks[a];
				m[Dy][Dx] ^= tanks[a];
				t[a] = (Dx << 4) | Dy;
				return false;
			}
			return false;
		}
		bool ugoku(int r[4],int ply) {
			u8 new_dick = dick & 0xf0;
			u8 om[9][9];
			memcpy(om, m, sizeof m);
			for (int i = ply<<1; i <(ply<<1)+2; ++i)
				if (!((dick >> i + 4) & 1))
					if (move_move(r[i])) {
						if (do_move(i, r[i], om))
							return true;
					}
			memcpy(om, m, sizeof m);
			for (int i = 0; i < 4; ++i)
				if (!((dick >> i) & 1) && !move_move(r[i])) {
					
					/*
					int x = t[i] >> 4, y = t[i] & 15;
					int tx = dx[r[i] - 4], ty = dy[r[i] - 4];
					int rz = om[y][x];
					*/
					new_dick |= 1u << i;
					/*
					for (; in_range(x += tx) && in_range(y += ty);)
						if (int c = om[y][x] & stopping) {
							switch (om[y][x] & 7) {
								case 0:
									//	cout<<"Tak"<<i<<" "<<mono(c)<<" "<<mono(rz)<<endl;
									if (mono(c) && mono(rz)) {
										int tid = TankId(c);
										if ((r[tid] ^ r[i]) != 2) {
											m[y][x] = 0, new_dick |= c << 1;
											//		cout<<"Fa?"<<endl;
										} else {}
										//		cout<<"Kie"<<endl;
									} else {
										m[y][x] = 0;
										new_dick |= c << 1;
									}
									break;
								case brick: m[y][x] = 0;
									break;
								case steel: break;
								case headq: return true;
								default: exit(1);
							}
							break;
						}
					*/
				}
			dick = new_dick | (new_dick >> 4);
			return false;
		};
		struct Map{
			int x[3];
			void clear(){x[0]=x[1]=x[2]=0;}
			void set(int X,int Y){
				x[0]=x[1]=x[2]=0;
				x[Y/3]|=1<<(Y%3*9+X);
			}
			bool get(int X,int Y){
				return (x[Y/3]>>(Y%3*9+X))&1;
			} 
			//((i * 27 + j) % 9, (i * 27 + j) / 9)
			//y*9+x)/27 (y*9+x)%27
			Map operator |(Map k){
				Map re;
				re.x[0]=x[0]|k.x[0];
				re.x[1]=x[1]|k.x[1];
				re.x[2]=x[2]|k.x[2];
				return re;
			}		
			Map operator &(Map k){
				Map re;
				re.x[0]=x[0]&k.x[0];
				re.x[1]=x[1]&k.x[1];
				re.x[2]=x[2]&k.x[2];
				return re;
			}
			Map move(int p){
				Map re;
				if(p==0){
					re.x[0]=x[0]>>9;
					re.x[0]|=(x[1]&511)<<18;
					re.x[1]=x[1]>>9;
					re.x[1]|=(x[2]&511)<<18;
					re.x[2]=x[2]>>9;
				}
				if(p==1){
					re.x[0]=(x[0]&133955070)>>1;
					re.x[1]=(x[1]&133955070)>>1;
					re.x[2]=(x[2]&133955070)>>1;
				}
				if(p==2){
					re.x[0]=(x[0]&262143)<<9;
					re.x[1]=(x[1]&262143)<<9;
					re.x[1]|=(x[0]&133955584)>>18;
					re.x[2]=(x[2]&262143)<<9;
					re.x[2]|=(x[1]&133955584)>>18;
				}
				if(p==3){
					re.x[0]=(x[0]&66977535)<<1;
					re.x[1]=(x[1]&66977535)<<1;
					re.x[2]=(x[2]&66977535)<<1;
				}
				return re;
			}
			void reverse(){
				int cnt=0,a[27];
				for(int I=0;I<3;I++){
					for(int i=0;i<27;i++)a[i]=x[I]&1,x[I]>>=1;
					for(int i=0;i<27;i++)x[I]<<=1,x[I]|=a[i];
				}
				swap(x[0],x[2]);
			}
			operator bool(){return x[0]|x[1]|x[2];}
		}fore,possible_position[2],x[9][9];
		int init(istream &z) {
			Json::Value r;
			Json::FastWriter fw;
			z >> r;
			memset(m, 0, sizeof m);
			m[0][4] = headq;
			m[8][4] = headq;
			memcpy(t, inilis, 4u);
			for (int i = 0; i < 4; ++i)m[t[i] & 15][t[i] >> 4] = tanks[i];

			auto c = r["requests"], d = r["responses"];
			auto init = c[0u];
			int ply = init["mySide"].asInt();
			if(ply){
				possible_position[0].set(t[0]>>4,t[0]&15);
				possible_position[1].set(t[1]>>4,t[1]&15);
			}else{
				possible_position[0].set(t[2]>>4,t[2]&15);
				possible_position[1].set(t[3]>>4,t[3]&15);
			}
			for(int i=0;i<9;i++)
				for(int j=0;j<9;j++){
					x[i][j].x[i/3]|=511<<(i%3*9);
					x[j][i].x[0]|=262657<<i;
					x[j][i].x[1]|=262657<<i;
					x[j][i].x[2]|=262657<<i;
				}
			auto ms = [this, &fw](Json::Value fuck, u8 shit) {
				for (int i = 0; i < 3; ++i) {
					int v = fuck[i].asInt();
					for (int j = 0; j < 27; ++j)
						if ((v >> j) & 1)m[i * 3 + j / 9][j % 9] = shit;
				}// i*3+j/9 j%9
				//  y/3 x+y%3*9
			};
			ms(init["brickfield"], brick);
			ms(init["steelfield"], steel);
			ms(init["waterfield"], water);
			for(int i=0;i<3;i++)fore.x[i]=init["forestfield"][i].asInt();
			int e = d.size();
			dick = 0;
			for (int i = 0; i < e; ++i) {
				
				//	cout<<"Step "<<i<<endl;
				//	Output();
				int op[4];
				Json::Value v1 = c[i+1]["action"],v2 = d[i];
				if (!ply)swap(v1, v2);
				op[0]=v1[0u].asInt();
				op[1] = v1[1u].asInt();
				op[2] = v2[0u].asInt(), op[3] = v2[1u].asInt();
				if (!ply)swap(v1, v2);
				v1=c[i+1];
				for(int j=0;j<v1["destroyed_blocks"].size();j+=2){
					int x = v1["destroyed_blocks"][j].asInt();
					int y = v1["destroyed_blocks"][j+1].asInt();
					m[y][x]=0;
				}
				for(int j=0;j<v1["destroyed_tanks"].size();j+=2){
					int x = v1["destroyed_tanks"][j].asInt();
					int y = v1["destroyed_tanks"][j+1].asInt();
					int c=m[y][x]&0x78;
					if(c==0)continue;
					int tid=TankId(c); 
					m[y][x] = 0, dick |=(c<<1)||(c>>3);
				}		
				if (ugoku(op,ply)) {
					cerr << "data parsing error" << endl;
					exit(1);
				}
				for(int j=0;j<v1["final_enemy_positions"].size()/2;j++){
					int x = v1["final_enemy_positions"][j<<1].asInt();
					int y = v1["final_enemy_positions"][(j<<1)|1].asInt();
					if(!ply)j+=2;
					if((x>=0)&&(dick&(1<<(j+4)))){
						dick&=255^(1<<(j|4))^(1<<j),m[y][x]^=tanks[j];
						t[j]=(x<<4)|y;
					}else if((x<0)&&((dick&(1<<(j+4)))==0)){
						int ty=t[j]&15,tx=t[j]>>4;
						dick|=(1<<(j+4))^(1<<j),m[ty][tx]^=tanks[j];
					}else if(x>=0){
						int ty=t[j]&15,tx=t[j]>>4;
						m[ty][tx]^=tanks[j];
						m[y][x]^=tanks[j];
						t[j]=(x<<4)|y;
					}
					if(!ply)j-=2;
//					printf("%d %d\n",x,y);
					if(x>=0)possible_position[j].set(x,y);
					else if(x==-1)possible_position[j].clear();
					else{
						possible_position[j]=possible_position[j]|possible_position[j].move(0)|possible_position[j].move(1)|possible_position[j].move(2)|possible_position[j].move(3);
						possible_position[j]=possible_position[j]&fore;
					}
				}
				//	for(int i=0;i<4;++i)cout<<out[op[i]+1]<<" ";cout<<endl;
			}
			/*
			for(int i=0;i<9;i++){
				for(int j=0;j<9;j++)printf("%d",possible_position[0].get(j,i));
				puts("");
			}
			puts("");
			for(int i=0;i<9;i++){
				for(int j=0;j<9;j++)printf("%d",possible_position[1].get(j,i));
				puts("");
			}
			puts("");
			Output();*/
			return ply;
		}
		inline bool canMove(int tankId) { return !((dick >> tankId + 4) & 1); }
		inline bool canShot(int tankId) { return !((dick >> tankId) & 1); }
		inline int tankX(int tankId) { return t[tankId] >> 4; }
		inline int tankY(int tankId) { return t[tankId] & 15; }
		void Output() {
			char DrawBoard[30][30];
			for (int i = 0; i < 30; ++i)
				for (int j = 0; j < 30; ++j) {
					char x = (i % 3 == 0) + (j % 3 == 0) * 2;
					switch (x) {
						case 0: x = ' ';
							break;
						case 1: x = '-';
							break;
						case 2: x = '|';
							break;
						case 3: x = '+';
							break;
					}
					DrawBoard[i][j] = x;
				}
			for (int i = 0; i < 4; ++i) {
				char x = (t[i] >> 4) + '0';
				char y = (t[i] & 15) + 'A';
				printf("tank %c(%c%c):%s%s\n", i + 'a', x, y, canMove(i) ? " move" : "", canShot(i) ? " shoot" : "");
			}
			for (int i = 0; i < 9; ++i) {
				DrawBoard[0][i * 3 + 1] = '0' + i;
				DrawBoard[i * 3 + 1][0] = 'A' + i;
			}
			for (int i = 0; i < 9; ++i)
				for (int j = 0; j < 9; ++j) {
#define __(a, b) DrawBoard[i*3+a][j*3+b]
					switch (m[i][j] & 7) {
						case space: // space or tank
							__(1, 1) = __(1, 2) = __(2, 1) = __(2, 2) = ' ';
							for (int z = 0; z < 4; ++z) {
								int x = 1 + (z >> 1), y = 1 + (z & 1);
								if (m[i][j] & tanks[z])
									__(x, y) = 'a' + z;
							}
							break;
						case water: __(1, 1) = __(1, 2) = __(2, 1) = __(2, 2) = '~';
							break;
						case brick: __(1, 1) = __(1, 2) = __(2, 1) = __(2, 2) = '.';
							break;
						case steel: __(1, 1) = __(1, 2) = __(2, 1) = __(2, 2) = '*';
							break;
						case headq: __(1, 1) = __(1, 2) = __(2, 1) = __(2, 2) = '$';
							break;
					}
				}
#undef __
			for (int i = 0; i < 28; ++i) {
				DrawBoard[i][28] = 0;
				puts(DrawBoard[i]);
			}
		}
		void swapSide() {
			for (int i = 0; i < 81; ++i) {
				int x = i % 9, y = i / 9;
				m[y][x] = (m[y][x] & 7) | ((m[y][x] & 24) << 2) | ((m[y][x] & 96) >> 2);
			}
			for (int i = 0; i < 41; ++i) {
				int x = i % 9, y = i / 9;
				swap(m[y][x], m[8 - y][8 - x]);
			}
			possible_position[0].reverse();
			possible_position[1].reverse();
			swap(t[0], t[2]);
			swap(t[1], t[3]);
			for (int i = 0; i < 4; ++i)
				t[i] = 136 - t[i];
			dick = ((dick & 51) << 2) | ((dick & 204u) >> 2);
		}
	};
	struct DistA {
		u8 mat[81], que[300];
		u8 stt;
		pos sttt;
		board *ref;
		void operator()(board &gg, u8 start) {
			
			const u8 ccmask = start == 4 ? tanks[2] | tanks[3] : (tanks[0] | tanks[1]);
			const u8 cdmask = start == 4 ? tanks[0] | tanks[1] : (tanks[2] | tanks[3]);
			ref = &gg;
			const bool Ending = ([&]() {
				if (start == 4)return gg.ise(2) || gg.ise(3);
				return gg.ise(0) || gg.ise(1);
			})();
			const int thisTank = start == 4 ? 2 : 0;
			memset(mat, 0xff, sizeof mat);
			int ql = 0, qr = 0;
			stt = que[qr++] = start;
			mat[start] = 0;
			int x = start % 9, y = start / 9;
			sttt = makepos(x, y);
			for (int i = 0; i < 4; ++i) {
				for (int tx = x, ty = y, cz = 1; in_range(tx += dx[i]) && in_range(ty += dy[i]);) {
					if (gg(tx, ty) == steel) break;
					if (gg(tx, ty) == headq)break;
					if(Ending && (gg(tx,ty)&cdmask) && (start==4?ty<3:ty>5))break;
					mat[que[qr++] = ty * 9 + tx] = cz, cz += (gg(tx, ty) == brick ? 2 : 0);
				}
			}
			while (ql != qr) {
				int c = que[ql++];
				int x = c % 9, y = c / 9, ds = mat[c];
				if (gg(x, y) & ccmask)continue;
				if (Ending)
					if (gg(x, y) & cdmask)
						if (start == 4) { if (y < 3)continue; }
						else if (y > 5)continue;
				for (int i = 0; i < 4; ++i) {
					int tx = x + dx[i], ty = y + dy[i];
					if ((!in_range(tx)) || (!in_range(ty)))continue;
					if (gg(tx, ty) == headq)continue;
					if (exdis[gg.m[ty][tx] & 7] == 255)continue;
					int tds = ds + 1 + (exdis[gg.m[y][x] & 7]);// + (gg.m[y][x] & tankmask ? 1 : 0);
					int tc = tx + ty * 9;
					if (tds < (int) mat[tc])mat[tc] = tds, que[qr++] = tc;
				}
			}
			if(start==4){			
				int u=-1;
				for(int i=0;i<81;i++)if(gg.possible_position[0].get(i%9,i/9)){
					if(u==-1)u=i;
					else if(mat[i]<mat[u])u=i;
				}
				if(u!=-1){
					int x=u%9,y=u/9;
					gg.t[2]=(x<<4)|y;
					gg.m[y][x]|=tanks[2];
					gg.dick&=255^4^64;
				}
				u=-1;
				for(int i=0;i<81;i++)if(gg.possible_position[1].get(i%9,i/9)){
					if(u==-1)u=i;
					else if(mat[i]<mat[u])u=i;
				}
				if(u!=-1){
					int x=u%9,y=u/9;
					gg.t[3]=(x<<4)|y;
					gg.m[y][x]|=tanks[3];
					gg.dick&=255^8^128;
				}
			//	gg.Output(); 
			}
		}
		void operator()(board &gg, pos t) { operator()(gg, t.x + t.y * 9); }
		void print() {
			for (int y = 0; y < 9; ++y)
				for (int x = 0; x < 9; ++x)
					printf("%03d%c", mat[y * 9 + x], x == 8 ? '\n' : ' ');
		}
		u8 operator()(int x, int y) {
			if (!in_range(x))return 255u;
			if (!in_range(y))return 255u;
			return mat[y * 9 + x];
		}
		u8 operator[](pos t) { return mat[t.y * 9 + t.x]; }
		bool inAtkRange(pos t) {
			if (t.x != sttt.x && t.y != sttt.y)return 0;
			return !ref->containSteel(t, sttt);
		}
		u8 shootBoom(pos t, int st) {
			pos r = ref->shopo(t, st);
			if (!valid(r))return 0;
			if ((ref->operator[](r) & 7) == steel)return 0;
			if ((ref->operator[](r) & 7) == 4 && operator[](r))return 0;
			if (inAtkRange(t) && inAtkRange(r) && (operator[](r) < operator[](t)))return 1;
			if (operator[](r) + maha(t, r) >= operator[](t))return 0;
			for (; !(t == r); t += dire[st])
				if (operator[](t + dire[st]) >= operator[](t))
					return 0;
			return 1;
		}
		u8 approx(pos t, int mov) {
			// return the approximate distance after taking a move
			if (mov == -1)return operator[](t);
			if (mov < 4)return operator[](movedir(t, mov));
			return operator[](t) - shootBoom(t, mov - 4);
		}
	};
	struct DistB {
		u8 mat[81], que[300];
		u8 stt;
		pos sttt;
		board *ref;
		void operator()(board &gg, const u8 start) {
			ref = &gg;
			memset(mat, 0xff, sizeof mat);
			int ql = 0, qr = 0;
			stt = que[qr++] = start;
			mat[start] = 0;
			while (ql != qr) {
				int c = que[ql++];
				int x = c % 9, y = c / 9, ds = mat[c];
				for (int i = 0; i < 4; ++i) {
					int tx = x + dx[i], ty = y + dy[i];
					if (!valid(makepos(tx, ty)))continue;
					if (gg(tx, ty) == headq)continue;
					if (exdis[gg.m[ty][tx] & 7] == 255)continue;
					int tds = ds + 1 + (exdis[gg.m[y][x] & 7]);// + (gg.m[y][x] & tankmask ? 1 : 0);
					int tc = tx + ty * 9;
					if (tds < (int) mat[tc])mat[tc] = tds, que[qr++] = tc;
				}
			}
		}
		void operator()(board &gg, pos t) { operator()(gg, t.x + t.y * 9); }
		void print() {
			for (int y = 0; y < 9; ++y)
				for (int x = 0; x < 9; ++x)
					printf("%03d%c", mat[y * 9 + x], x == 8 ? '\n' : ' ');
		}
		u8 operator()(int x, int y) {
			if (!in_range(x))return 255u;
			if (!in_range(y))return 255u;
			return mat[y * 9 + x];
		}
		u8 operator[](pos t) { return mat[t.y * 9 + t.x]; }
	};
}
#include "string.h"
#include "stdint.h"
#include <cassert>
namespace supp {
	using u8=uint8_t;
	template<class T, u8 siz>
	struct svec {
		static const u8 size = siz;
		using content=T;
		using this_type=svec<content, size>;
		T buf[siz];
		u8 len;
		svec() : len(0) { memset(buf, 0, sizeof buf); }
		svec(T *a, T *b) : len(b - a) {
			assert(b - a <= siz);
			memcpy(buf, a, (b - a) * sizeof(T));
		}
		template<class Z>
		T &operator[](Z t) {
			assert(t >= 0 && t < len);
			return buf[t];
		}
		this_type &push(T a) {
			assert(len < siz);
			buf[len++] = a;
			return *this;
		}
		template<class fn>
		this_type filter(fn c) const {
			this_type ne;
			for (int i = 0; i < len; ++i)
				if (c(buf[i]))
					ne.push(buf[i]);
			return ne;
		}
		template<class fn>
		auto map(fn c) const {
			svec<decltype(c(buf[0])), siz> r;
			for (int i = 0; i < len; ++i)
				r.push(c(buf[i]));
			return r;
		}
		template<class fn, class ty>
		ty reduce(fn c, ty r) const {
			for (int i = 0; i < len; ++i)
				r = c(r, buf[i]);
			return r;
		}
		template<class fn>
		this_type minvec(fn func) const {
			if (!len)return this_type();
			this_type ne;
			auto r = func(buf[0]);
			ne.push(buf[0]);
			for (int i = 1; i < len; ++i) {
				auto c = func(buf[i]);
				if (c < r) {
					ne.len = 0;
					r = c;
				}
				if (c == r) { ne.push(buf[i]); }
			}
			return ne;
		}
		operator bool() const { return len; }
	};
}
#include "cassert"
#define _log(msg) cerr<<msg<<endl;
namespace Policy {
	using namespace Complete_Simulator;
	struct Walk {
		using i8=int8_t;
		using opvec=supp::svec<i8, 9>;
		using dist=Complete_Simulator::DistA;
		dist myh, opph;
		Complete_Simulator::board *glo;
		Walk(Complete_Simulator::board *gl) : glo(gl) {
			myh(*glo, 4);
			opph(*glo, 8 * 9 + 4);
		}
		opvec allMove(int tank, int letShoot = 1, int letStill = 1) {
			const pos t = glo->tank(tank);
			const bool canS = letShoot && glo->canShot(tank);
			opvec ret;
			if (letStill) ret.push(-1);
			for (int i = 0; i < 4; ++i)
				if (valid(t + dire[i])) {
					if (!(glo->operator[](movedir(t, i))))ret.push(i);
					if (canS) ret.push(i + 4);
				}
			return ret;
		}
		opvec allSafeMove(int tank, int letShoot = 1) {
			const pos thi = glo->tank(tank);
			return allMove(tank, letShoot).filter([this, tank, thi](int move) {
				return glo->dangerous(tank, move) == 0;
			});
		}
		opvec filterAttack(opvec e, int tank, dist &c) {
			const pos thi = glo->tank(tank);
			const int dd = opph[thi];
			return e.filter([&](int move) {
				//	cerr << (int) move << "[:]" << dd << "->" << (int) opph.approx(thi, move) << " ";
				return opph.approx(thi, move) < dd;
			});
		}
		i8 randOp(opvec c) {
			const int len = c.len;
			if (len == 0) return -1;
			return c[unsigned_sample(rnpool, len)];
		}
		i8 shootOpp(int tank) {
			const pos thi = glo->tank(tank);
			for (int i = 0; i < 4; ++i) {
				const pos cc = glo->shopo(thi, i);
				//	cerr << (int) cc.x << ">>" << (int) cc.y << endl;
				if ((*glo)[cc] & (tanks[2] | tanks[3]))
					return i + 4;
			}
			return 0;
		}
		i8 filterMinimum(int tank, opvec a) {
			const pos thi = glo->tank(tank);
			//		a.map([&](i8 mov){cout<<(int)mov<<":"<<(int)opph.approx(thi,mov)<<" ";return 0;});
			return randOp(a.minvec([thi, this](int mov) { return int(opph.approx(thi, mov) * 2) + (mov > 3); }));
		}
		i8 operateSemi(int tank, int move) {
			//	cerr<<tank<<" semi "<<move<<endl;
			const bool canS = glo->canShot(tank);
			const pos thi = glo->tank(tank);
			const int s2 = glo->dangerous(tank, move, tank ^ 2) == 1 ? glo->TypeOfSemi(tank, move, tank ^ 2) : 0;
			const int s3 = glo->dangerous(tank, move, tank ^ 3) == 1 ? glo->TypeOfSemi(tank, move, tank ^ 3) : 0;
			/** Need to decide type of semidan.
			/// 1. XX Me
			///    Op
			///    Shoot Move Stop
			///    14    3    6       ----> 0 1 7
			///          1    3       ----> 0 1 7
			/// 2. XX XX Me
			///    Op
			///    Move Stop
			 *    1    7
			 *
			 * 3. XX    XX
			 *    XX MO XX
			 *    XX    XX
			 *    A: Attack
			 *    G: Guard
			 *    M: Move
			 *    S: Shoot
			 *    S  Stop
			 *    AM AS S GM GS
			 *    5       3  1  (Same dis)
			 *    4       3
			 *    2 0  3  0  0  (Near)
			 *    2 0  3  0  0
			 *    2 2  0  2  1  (Far)
			 *    0       1
			 */
			const int mainPivot = s2 ? s3 ? min(s2, s3) : s2 : s3;
			const int dire = move > 3 ? move - 4 : move;
			const int my = opph[thi], op = myh[thi];
			const int cr = (my <= op ? my == op ? 0 : 2 : 4) + (canS);
			const opvec directions = allMove(tank, 0, 0);
			const int AD = randOp(directions.minvec([&](int mov) { return opph[movedir(thi, mov)]; }));
			const int GD = randOp(directions.minvec([&](int mov) { return myh[movedir(thi, mov)]; }));
			//	cerr<<s2<<" "<<s3<<"=>"<<mainPivot<<endl;
			//	cerr<<"AD "<<AD<<" GD "<<GD<<" cr "<<cr<<endl;
			array<int, 4> cr1, cr2;
			switch (mainPivot) {
				case 1: if (canS)return select(3, cr1 = {0, 1, 7}, cr2 = {dire + 4, dire, -1});
					return select(2, cr1 = {1, 7}, cr2 = {dire, -1});
				case 2: return select(2, cr1 = {1, 7}, cr2 = {dire, -1});
				case 3:
					switch (cr) {
						case 0: return select(2, cr1 = {4, 3}, cr2 = {AD, GD});
						case 1: return select(3, cr1 = {3, 3, 1}, cr2 = {AD, GD, GD + 4});
						case 2:
						case 3: return select(2, cr1 = {2, 3}, cr2 = {AD, -1});
						case 4: return GD;
						case 5: return select(4, cr1 = {2, 2, 2, 1}, cr2 = {AD, AD + 4, GD, GD + 4});
						default: assert(0);
					}
					break;
				default:assert(0);
			}
		}
		i8 genMove(int tank) {
			if (!glo->canMove(tank))return -1;
			const pos thi = glo->tank(tank),
				tha = glo->tank(tank ^ 1),
				t1 = glo->tank(tank ^ 2),
				t2 = glo->tank(tank ^ 3);
			if (glo->canShot(tank) && opph[thi] == 1)
				return (glo->shoodire(thi, makepos(4, 8)) ^ 2) + 4;
			const int is_left = thi.x < tha.x ? 1 : (thi.x == tha.x ? tank & 1 : 0);
			const int ris_left = t1.x < t2.x ? 1 : (t1.x == t2.x ? tank & 1 : 0);
			const int oppCh =
				abs(t1.x - (int) thi.x) < abs(t2.x - (int) thi.x) ? 0 : (abs(t1.x - (int) thi.x) == abs(t2.x - (int) thi.x)
																							? (is_left == ris_left ? 0 : 1) : 1);
			int oppTnk = oppCh ? (tank ^ 3) : (tank ^ 2);
			if (!glo->canMove(oppTnk))oppTnk ^= 1;
			const pos oppT = glo->tank(oppTnk);
			const i8 mov1 = thi.x <= 4 ? 1 : 3, sho1 = thi.x <= 4 ? 5 : 7;
			_log("Tank " << tank)
			if (opph[thi] == 1) {
				const int t = glo->shotTo(thi, home[1]);
				if ((!glo->dangerous(tank, -1)) || glo->dangerous(tank, t))return -1;
				const pos tt = movedir(thi, t);
				if ((*glo)[tt])return -1;
				return t;
			}
			/// Sun dog biss
			
			if (!glo->canShot(tank ^ 2) && !allMove(tank ^ 2).filter([&](int mov) {
				return glo->dangerous(tank ^ 2, mov, tank) != 2;
			}))
				return glo->shotTo(thi, t1) + 4;
			if (!glo->canShot(tank ^ 3) && !allMove(tank ^ 3).filter([&](int mov) {
				return glo->dangerous(tank ^ 3, mov, tank) != 2;
			}))
				return glo->shotTo(thi, t2) + 4;
			
			/// Gang!
			if (thi.y < 4 && oppT.y > 4 && opph[thi] > myh[oppT]) {
				const opvec cr = allSafeMove(tank);
				DistB disk;
				disk(*glo, oppT);
				//	cerr<<endl;
				//	disk.print();
				return randOp(cr.minvec([&](int mov) {
					if (mov < 4) { return int(disk[movedir(thi, mov)])*2; }
					if (glo->operator[](movedir(thi, mov - 4)) == brick)
						return (int(disk[movedir(thi, mov - 4)]) + 1)*2+1;
					else
						return (int(disk[thi]))*2+1;
				}));
			}
			_log("_")
			if (thi.y == 8 && oppT.y == 8 && glo->dangerous(tank, -1)) return 0;
			if (true) {
				_log("Ba")
			_l_attack:
				_log("La")
				const opvec safe = allSafeMove(tank);
				const opvec safeA = filterAttack(safe, tank, opph);
				safe.map([&](int mov) {
		//			cerr << mov << "(:)" << (int) opph[thi] << "," << (int) opph.approx(thi, mov) << " ";
					return 0;
				});
				_log("Lb")
				if (safeA.len)return filterMinimum(tank, safeA);
				_log("Lc")
				if (const opvec hvec = allMove(tank, 0).filter([&](int mov) {
					//			cerr << (int) mov << ":" << glo->dangerous(tank, mov, tank ^ 2) << ","
					//				  << glo->dangerous(tank, mov, tank ^ 3) << " ";
					return 1 == glo->dangerous(tank, mov);
				})) {
					return operateSemi(tank, filterMinimum(tank, hvec));
				}
				if (thi.y >= oppT.y) {
					_log("La3")
					const opvec smv = safe.filter([thi, tank, this](int mov) {
						if (mov > 3)return false;
						if (opph[thi] == opph[movedir(thi, mov)])return true;
						return false;
					});
					if (smv.len)return filterMinimum(tank, smv);
					if (glo->canShot(tank)) {
						if (int t = shootOpp(tank))return t;
						_log("Laa")
						assert(0);
						return rnpool->extrude_bit(myh[oppT] - opph[thi]) == 0 ? filterMinimum(tank, allMove(tank))
																								 : shootOpp(tank);
					}
					_log("La4")
					if (safe.len)return filterMinimum(tank, safe);
					_log("La5")
					return filterMinimum(tank, allMove(tank));
				}
				_log("Ld")
				if (glo->canShot(tank))
					if (int t = shootOpp(tank))return t;
				_log("Le")
				if (safe)return filterMinimum(tank, safe);
				_log("Lf")
				return filterMinimum(tank, allMove(tank));
			}
		}
	};
}
using i8=signed char;
Complete_Simulator::board global;
int ply;
void read() {
	ply = global.init(cin);
	if (ply)global.swapSide();
}
int revert(int t) { return t >= 0 ? t ^ 2 : t; }
int main(int argc, const char **argv) {
//	freopen("in.txt","r",stdin);
//	freopen("out.txt","w",stdout);
	srandit();
	read();
//	global.Output();
	Policy::Walk c(&global);
//	c.opph.print();
	int a = c.genMove(0);
	int b = c.genMove(1);
	if (a > 3) {
		global[global.tank(1)] ^= Complete_Simulator::tanks[1];
		if (global.shopo(Complete_Simulator::movedir(global.tank(1), b), (a - 4) ^ 2) == global.tank(0))
			a = -1;
		global[global.tank(1)] ^= Complete_Simulator::tanks[1];
	}
	if (b > 3) {
		global[global.tank(0)] ^= Complete_Simulator::tanks[0];
		if (global.shopo(Complete_Simulator::movedir(global.tank(0), a), (b - 4) ^ 2) == global.tank(1))
			b = -1;
		global[global.tank(0)] ^= Complete_Simulator::tanks[0];
	}
	if (ply)
		cout << "[" << revert(a) << "," << revert(b) << "]" << endl;
	else
		cout << "[" << a << "," << b << "]" << endl;
	return 0;
}
