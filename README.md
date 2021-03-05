# Matrix-characteristic-value

使用幂法、QR法、雅克比法求矩阵的特征值

幂法求矩阵特征值说明： 幂法求特征值的每次迭代时，将向量归一化，迭代公式为 yk=xk/max{xk} xk+1=Ayk max{xk}为向量 xk的无穷范数（最大分量）

Jacobi 方法说明： 首先初始化特征向量为一个单位矩阵，主对角线的元素为 1，其余元素为 0. 在要求解特征值的矩阵中的非主对角线元素中，找到绝对值最大元素 Aij tan2θ=2Aij/(Ajj-Aii)，得出 cosθ，sinθ 令矩阵 A’与 A 是对称相似的，有： Aii’=(cosθ)^2Aii-2cosθsinθ Aij+(sinθ)^2Ajj Ajj’=(cosθ)^2Ajj+2cosθsinθ Aij+(sinθ)^2Aii Aij’=Aji’=(cosθ-sinθ)Aij+cosθsinθ(Aii-Ajj) 若当前迭代的矩阵的非主对角线元素中最大值小于给定的精度时，停止迭代, 否则令 A =A’, 进行迭代，最终得到特征值以及特征向量

QR 方法说明： 首先使用高斯消元法求得矩阵的上 Hensenberg 矩阵 对求得的上 Hensenberg 矩阵进行 QR 分解，QR 分解的迭代公式为： Ak=QkRk Ak+1=Rk+1Qk+1
多次计算 A，直到 A 的变化小于给定的精度时，停止迭代，此时 A 主对角线上的元素即为矩 阵的全部特征值
