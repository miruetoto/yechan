---
title: (책공부) [SOGT] General Topology
author: 신록예찬
date: 08/31/2025
draft: false
---

## Chap 5. Topological Spaces: Definitions

### TOPOLOGICAL SPACES

$X$를 공집합이 아닌 집합이라 하자. $X$의 부분집합들의 모임 $\mathcal{T}$가 $X$ 위의 위상이란, $\mathcal{T}$가 다음 공리들을 만족하는 것이다.

- [$O_1$] $X$ and $\varnothing$ belong to $\mathcal{T}$.
- [$O_2$] The union of any number of sets in $\mathcal{T}$ belongs to $\mathcal{T}$.
- [$O_3$] The intersection of any two sets in $\mathcal{T}$ belongs to $\mathcal{T}$.

$\mathcal{T}$의 원소들을 $\mathcal{T}$-열린 집합, 또는 간단히 열린 집합이라 부르며, $X$와 $\mathcal{T}$를 함께 쌍 $(X, \mathcal{T})$를 위상 공간이라 한다.

Example 1.1: $\mathcal{U}$를 4장에서 다룬 실수의 모든 열린 집합들의 모임이라 하자. 그러면 $\mathcal{U}$는 $\mathbb{R}$ 위의 위상이며, 이를 $\mathbb{R}$ 위의 보통위상이라 한다. 마찬가지로 평면 $\mathbb{R}^2$에서 모든 열린 집합들의 모임 $\mathcal{U}$도 위상이며, $\mathbb{R}^2$ 위의 보통위상이라 한다. 별도로 명시하지 않는 한 $\mathbb{R}$과 $\mathbb{R}^2$에서는 항상 보통위상을 가정한다.

Example 1.2: $X = \{a, b, c, d, e\}$의 부분집합들의 다음 모임들을 생각하자.
$$\mathcal{T}_1 = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{b, c, d, e\}\}$$
$$\mathcal{T}_2 = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{b, c, d\}\}$$
$$\mathcal{T}_3 = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{a, b, d, e\}\}$$

$\mathcal{T}_1$은 세 공리 [$O_1$], [$O_2$], [$O_3$]을 만족하므로 $X$ 위의 위상임을 확인하라. 그러나 $\mathcal{T}_2$는 $X$ 위의 위상이 아닌데, $\mathcal{T}_2$의 두 원소의 합집합
$$\{a, c, d\} \cup \{b, c, d\} = \{a, b, c, d\}$$
이 $\mathcal{T}_2$에 속하지 않기 때문이다. 즉 $\mathcal{T}_2$는 공리 [$O_2$]를 만족하지 않는다.

또한 $\mathcal{T}_3$도 $X$ 위의 위상이 아닌데, $\mathcal{T}_3$의 두 집합의 교집합
$$\{a, c, d\} \cap \{a, b, d, e\} = \{a, d\}$$
이 $\mathcal{T}_3$에 속하지 않기 때문이다. 즉 $\mathcal{T}_3$은 공리 [$O_3$]을 만족하지 않는다.

Example 1.3: $\mathcal{D}$를 $X$의 모든 부분집합들의 모임이라 하자. $\mathcal{D}$가 $X$ 위의 위상에 대한 공리들을 만족함을 확인하라. 이 위상을 이산위상이라 하며, $X$와 그 이산위상을 함께, 즉 쌍 $(X, \mathcal{D})$를 이산위상공간 또는 간단히 이산공간이라 한다.

Example 1.4: 공리 [$O_1$]에서 보듯이, $X$ 위의 위상은 반드시 집합 $X$와 $\varnothing$를 포함해야 한다. $X$와 $\varnothing$만으로 이루어진 모임 $\mathcal{J} = \{X, \varnothing\}$은 그 자체로 $X$ 위의 위상이다. 이를 비이산위상이라 하며, $X$와 그 비이산위상을 함께, 즉 $(X, \mathcal{J})$를 비이산위상공간 또는 간단히 비이산공간이라 한다.

Example 1.5: $\mathcal{T}$를 여집합이 유한인 $X$의 모든 부분집합들의 모임에 공집합 $\varnothing$를 함께 넣은 것이라 하자. 이 모임 $\mathcal{T}$도 $X$ 위의 위상이다. 이를 여유한위상 또는 $X$ 위의 $T_1$-위상이라 한다. ($T_1$의 의미는 이후 장에서 다룬다.)

Example 1.6: $X$ 위의 임의의 두 위상 $\mathcal{T}_1$과 $\mathcal{T}_2$의 교집합 $\mathcal{T}_1 \cap \mathcal{T}_2$도 $X$ 위의 위상이다. [$O_1$]에 의해 $X$와 $\varnothing$는 각각 $\mathcal{T}_1$과 $\mathcal{T}_2$ 모두에 속하므로, $X$와 $\varnothing$는 교집합 $\mathcal{T}_1 \cap \mathcal{T}_2$에 속한다. 즉 $\mathcal{T}_1 \cap \mathcal{T}_2$는 [$O_1$]을 만족한다. 또한 $G, H \in \mathcal{T}_1 \cap \mathcal{T}_2$이면, 특히 $G, H \in \mathcal{T}_1$이고 $G, H \in \mathcal{T}_2$이다. 그런데 $\mathcal{T}_1$과 $\mathcal{T}_2$가 위상이므로 $G \cap H \in \mathcal{T}_1$이고 $G \cap H \in \mathcal{T}_2$이다. 따라서
$$G \cap H \in \mathcal{T}_1 \cap \mathcal{T}_2$$
즉 $\mathcal{T}_1 \cap \mathcal{T}_2$는 [$O_3$]을 만족한다. 마찬가지로 $\mathcal{T}_1 \cap \mathcal{T}_2$는 [$O_2$]를 만족한다.

앞의 예제에서의 진술은 사실 위상들의 임의의 모임으로 일반화할 수 있다. 즉,

Theorem 5.1: Let $\{\mathcal{T}_i : i \in I\}$ be any collection of topologies on a set $X$. Then the intersection $\bigcap_i \mathcal{T}_i$ is also a topology on $X$.

마지막 예제에서, 위상들의 합집합은 위상이 되지 않을 수 있음을 보인다.

Example 1.7: 각 모임
$$\mathcal{T}_1 = \{X, \varnothing, \{a\}\} \quad \text{and} \quad \mathcal{T}_2 = \{X, \varnothing, \{b\}\}$$
은 $X = \{a, b, c\}$ 위의 위상이다. 그러나 합집합
$$\mathcal{T}_1 \cup \mathcal{T}_2 = \{X, \varnothing, \{a\}, \{b\}\}$$
은 $X$ 위의 위상이 아닌데, [$O_2$]를 위반하기 때문이다. 즉 $\{a\} \in \mathcal{T}_1 \cup \mathcal{T}_2$, $\{b\} \in \mathcal{T}_1 \cup \mathcal{T}_2$이지만 $\{a\} \cup \{b\} = \{a, b\}$는 $\mathcal{T}_1 \cup \mathcal{T}_2$에 속하지 않는다.

$G$가 점 $p \in X$를 포함하는 열린 집합이면, $G$를 $p$의 열린 근방이라 한다. 또한 $p$를 제외한 $G$, 즉 $G \setminus \{p\}$를 $p$의 제거된 열린 근방이라 한다.

Remark: 공리 [$O_1$], [$O_2$], [$O_3$]은 다음 두 공리와 동치이다:

- [$O_1^*$] The union of any number of sets in $\mathcal{T}$ belongs to $\mathcal{T}$.
- [$O_2^*$] The intersection of any finite number of sets in $\mathcal{T}$ belongs to $\mathcal{T}$.

[$O_1^*$]이 $\varnothing$가 $\mathcal{T}$에 속함을 함의하는데, 이는
$$\cup \{G \in \mathcal{T} : G \in \varnothing\} = \varnothing$$
즉 집합들의 공합집합은 공집합이기 때문이다. 또한 [$O_2^*$]는 $X$가 $\mathcal{T}$에 속함을 함의하는데, 이는
$$\cap \{G \in \mathcal{T} : G \in \varnothing\} = X$$
즉 $X$의 부분집합들의 공교집합은 $X$ 자체이기 때문이다.

### ACCUMULATION POINTS

$X$를 위상 공간이라 하자. 점 $p \in X$가 $X$의 부분집합 $A$의 집적점 또는 극한점(군점 또는 도래점이라고도 함)이란, $p$를 포함하는 모든 열린 집합 $G$가 $p$와 다른 $A$의 점을 포함하는 것이다. 즉,
$$G \text{ open}, \ p \in G \quad \text{implies} \quad (G \setminus \{p\}) \cap A \neq \varnothing$$
$A$의 집적점들의 집합을 $A'$로 나타내며, $A$의 도래집합이라 한다.

Example 2.1: 모임 $\mathcal{T} = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{b, c, d, e\}\}$은 $X = \{a, b, c, d, e\}$ 위의 위상을 정의한다. $X$의 부분집합 $A = \{a, b, c\}$를 생각하자. $b \in X$는 $A$의 극한점인데, $b$를 포함하는 열린 집합들이 $\{b, c, d, e\}$와 $X$이고, 각각 $b$와 다른 $A$의 점, 즉 $c$를 포함하기 때문이다. 반면에 점 $a \in X$는 $A$의 극한점이 아닌데, $a$를 포함하는 열린 집합 $\{a\}$가 $a$와 다른 $A$의 점을 포함하지 않기 때문이다. 마찬가지로 점 $d$와 $e$는 $A$의 극한점이고, 점 $c$는 $A$의 극한점이 아니다. 따라서 $A' = \{b, d, e\}$가 $A$의 도래집합이다.

Example 2.2: $X$를 비이산 공간이라 하자. 즉 $X$와 $\varnothing$만이 $X$의 유일한 열린 부분집합이다. 그러면 $X$는 임의의 점 $p \in X$를 포함하는 유일한 열린 집합이다. 따라서 $p$는 공집합 $\varnothing$와 $p$만으로 이루어진 집합 즉 한원소 집합 $\{p\}$를 제외한 $X$의 모든 부분집합의 집적점이다. 따라서 $X$의 임의의 부분집합 $A$의 도래집합 $A'$은 다음과 같다:
$$A' = \begin{cases} \varnothing & \text{if } A = \varnothing \\ \{p\}^c = X \setminus \{p\} & \text{if } A = \{p\} \\ X & \text{if } A \text{ contains two or more points} \end{cases}$$

직선 $\mathbb{R}$과 평면 $\mathbb{R}^2$의 보통위상에 대해, 위의 집적점의 정의는 4장에서 주어진 것과 동일함을 확인하라.

### CLOSED SETS

$X$를 위상 공간이라 하자. $X$의 부분집합 $A$가 닫힌 집합이란, 그 여집합 $A^c$가 열린 집합인 것이다.

Example 3.1: 모임 $\mathcal{T} = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{b, c, d, e\}\}$은 $X = \{a, b, c, d, e\}$ 위의 위상을 정의한다. $X$의 닫힌 부분집합들은
$$\varnothing, \ X, \ \{b, c, d, e\}, \ \{a, b, e\}, \ \{b, e\}, \ \{a\}$$
즉 $X$의 열린 부분집합들의 여집합이다. $\{b, c, d, e\}$와 같이 열리면서 동시에 닫힌 $X$의 부분집합이 있고, $\{a, b\}$와 같이 열리지도 않고 닫히지도 않은 $X$의 부분집합도 있음에 주목하라.

Example 3.2: $X$를 이산 위상 공간이라 하자. 즉 $X$의 모든 부분집합이 열려 있다. 그러면 $X$의 모든 부분집합은 그 여집합이 항상 열려 있으므로 닫혀 있기도 하다. 다시 말해 $X$의 모든 부분집합은 열리면서 동시에 닫혀 있다.

공간 $X$의 임의의 부분집합 $A$에 대해 $A^{cc} = A$임을 상기하라. 따라서

Proposition 5.2: In a topological space $X$, a subset $A$ of $X$ is open if and only if its complement is closed.

위상 공간의 공리 [$O_1$], [$O_2$], [$O_3$]과 드 모르간 법칙으로부터 다음을 얻는다

Theorem 5.3: Let $X$ be a topological space. Then the class of closed subsets of $X$ possesses the following properties:

- (i) $X$ and $\varnothing$ are closed sets.
- (ii) The intersection of any number of closed sets is closed.
- (iii) The union of any two closed sets is closed.

닫힌 집합은 그 극한점을 이용하여 다음과 같이 특성화할 수도 있다:

Theorem 5.4: A subset $A$ of a topological space $X$ is closed if and only if $A$ contains each of its accumulation points.

다시 말해, 집합 $A$가 닫혀 있을 필요충분조건은 $A$의 도래집합 $A'$이 $A$의 부분집합인 것, 즉 $A' \subset A$이다.

### CLOSURE OF A SET

$A$를 위상 공간 $X$의 부분집합이라 하자. $A$의 폐포란,
$$\bar{A} \quad \text{or} \quad A^-$$
로 나타내며, $A$를 포함하는 모든 닫힌 상위집합들의 교집합이다. 다시 말해 $\{F_i : i \in I\}$가 $A$를 포함하는 $X$의 모든 닫힌 부분집합들의 모임이면,
$$\bar{A} = \cap_i F_i.$$
$\bar{A}$는 닫힌 집합들의 교집합이므로 닫힌 집합임을 먼저 확인하라. 더 나아가 $\bar{A}$는 $A$를 포함하는 가장 작은 닫힌 상위집합이다. 즉 $F$가 $A$를 포함하는 닫힌 집합이면
$$A \subset \bar{A} \subset F$$
따라서 집합 $A$가 닫혀 있을 필요충분조건은 $A = \bar{A}$이다. 이 결과들을 형식적으로 진술하면:

Proposition 5.5: Let $\bar{A}$ be the closure of a set $A$. Then: (i) $\bar{A}$ is closed; (ii) if $F$ is a closed superset of $A$, then $A \subset \bar{A} \subset F$; and (iii) $A$ is closed iff $A = \bar{A}$.

Example 4.1: Example 3.1의 $X = \{a, b, c, d, e\}$ 위의 위상 $\mathcal{T}$를 생각하자. 닫힌 부분집합들은
$$\varnothing, \ X, \ \{b, c, d, e\}, \ \{a, b, e\}, \ \{b, e\}, \ \{a\}$$
따라서,
$$\overline{\{b\}} = \{b, e\}, \quad \overline{\{a, c\}} = X, \quad \overline{\{b, d\}} = \{b, c, d, e\}.$$
Example 4.2: $X$를 여유한 위상 공간이라 하자. 즉 유한 집합의 여집합과 $\varnothing$가 열린 집합이다. 그러면 닫힌 집합들은 정확히 $X$의 유한 부분집합들과 $X$ 자체이다. 따라서 $A \subset X$가 유한이면, $A$ 자체가 닫혀 있으므로 그 폐포 $\bar{A}$는 $A$이다. 반면에 $A \subset X$가 무한이면, $X$가 $A$를 포함하는 유일한 닫힌 상위집합이므로 $\bar{A}$는 $X$이다. 더 간결하게, 여유한 공간 $X$의 임의의 부분집합 $A$에 대해,
$$\bar{A} = \begin{cases} A & \text{if } A \text{ is finite} \\ X & \text{if } A \text{ is infinite} \end{cases}$$
집합의 폐포는 그 극한점을 이용하여 다음과 같이 완전히 기술할 수 있다:

Theorem 5.6: Let $A$ be a subset of a topological space $X$. Then the closure of $A$ is the union of $A$ and its set of accumulation points, i.e.,
$$\bar{A} = A \cup A'$$

점 $p \in X$가 $A \subset X$의 폐포점 또는 부착점이란, $p$가 $A$의 폐포에 속하는 것, 즉 $p \in \bar{A}$이다. 위 정리에 의하면, $p \in X$가 $A \subset X$의 폐포점일 필요충분조건은 $p \in A$이거나 $p$가 $A$의 극한점인 것이다.

Example 4.3: 유리수의 집합 $\mathbb{Q}$를 생각하자. 앞서 보았듯이 $\mathbb{R}$의 보통위상에서 모든 실수 $a \in \mathbb{R}$은 $\mathbb{Q}$의 극한점이다. 따라서 $\mathbb{Q}$의 폐포는 실수 전체의 집합, 즉 $\bar{\mathbb{Q}} = \mathbb{R}$이다.

위상 공간 $X$의 부분집합 $A$가 $B \subset X$에서 조밀하다란, $B$가 $A$의 폐포에 포함되는 것, 즉 $B \subset \bar{A}$이다. 특히 $A$가 $X$에서 조밀하다, 또는 $X$의 조밀한 부분집합이란, $\bar{A} = X$인 것이다.

Example 4.4: Example 4.1에서 다음을 확인하라.
$$\overline{\{a, c\}} = X \quad \text{and} \quad \overline{\{b, d\}} = \{b, c, d, e\}$$
여기서 $X = \{a, b, c, d, e\}$이다. 따라서 집합 $\{a, c\}$는 $X$의 조밀한 부분집합이지만 $\{b, d\}$는 그렇지 않다.

Example 4.5: Example 4.3에서 언급했듯이 $\bar{\mathbb{Q}} = \mathbb{R}$이다. 다시 말해 보통위상에서 유리수의 집합 $\mathbb{Q}$는 $\mathbb{R}$에서 조밀하다.

"폐포" 연산자는, $X$의 각 부분집합 $A$에 그 폐포 $\bar{A} \subset X$를 대응시키며, 아래 명제에 나오는 네 가지 성질을 만족하는데, 이를 쿠라토프스키 폐포 공리라 한다. 사실 이 공리들은 이후에 증명하겠지만 $X$ 위의 위상을 정의하는 데 사용될 수 있다.

Proposition 5.7: (i) $\bar{\varnothing} = \varnothing$; (ii) $A \subset \bar{A}$; (iii) $\overline{A \cup B} = \bar{A} \cup \bar{B}$; and (iv) $\bar{\bar{A}} = \bar{A}$.

### INTERIOR, EXTERIOR, BOUNDARY

$A$를 위상 공간 $X$의 부분집합이라 하자. 점 $p \in A$가 $A$의 내점이란, $p$가 $A$에 포함된 어떤 열린 집합 $G$에 속하는 것이다:
$$p \in G \subset A \quad \text{where } G \text{ is open}$$
$A$의 내점들의 집합을
$$\text{int}(A), \quad \mathring{A} \quad \text{or} \quad A^\circ$$
로 나타내며, $A$의 내부라 한다. $A$의 내부는 다음과 같이 특성화할 수도 있다:

Proposition 5.8: The interior of a set $A$ is the union of all open subsets of $A$. Furthermore: (i) $A^\circ$ is open; (ii) $A^\circ$ is the largest open subset of $A$, i.e. if $G$ is an open subset of $A$ then $G \subset A^\circ \subset A$; and (iii) $A$ is open iff $A = A^\circ$.

$A$의 외부란, $\text{ext}(A)$로 쓰며, $A$의 여집합의 내부, 즉 $\text{int}(A^c)$이다. $A$의 경계란, $b(A)$로 쓰며, $A$의 내부에도 외부에도 속하지 않는 점들의 집합이다. 다음으로 내부, 외부, 폐포 사이의 중요한 관계가 이어진다.

Theorem 5.9: Let $A$ be any subset of a topological space $X$. Then the closure of $A$ is the union of the interior and boundary of $A$, i.e. $\bar{A} = A^\circ \cup b(A)$.

Example 5.1: 끝점이 $a$와 $b$인 네 구간 $[a, b]$, $(a, b)$, $(a, b]$, $[a, b)$를 생각하자. 각각의 내부는 열린 구간 $(a, b)$이고, 각각의 경계는 끝점들의 집합, 즉 $\{a, b\}$이다.

Example 5.2: $X = \{a, b, c, d, e\}$ 위의 위상
$$\mathcal{T} = \{X, \varnothing, \{a\}, \{c, d\}, \{a, c, d\}, \{b, c, d, e\}\}$$
과 $X$의 부분집합 $A = \{b, c, d\}$를 생각하자. 점 $c$와 $d$는 각각 $A$의 내점인데,
$$c, d \in \{c, d\} \subset A$$
이고 $\{c, d\}$는 열린 집합이기 때문이다. 점 $b \in A$는 $A$의 내점이 아니다. 따라서 $\text{int}(A) = \{c, d\}$이다. $X$의 점 중 오직 $a$만이 $A$에 대해 외점, 즉 $A$의 여집합 $A^c = \{a, e\}$의 내점이다. 따라서 $\text{int}(A^c) = \{a\}$이다. 따라서 $A$의 경계는 점 $b$와 $e$로 이루어진다. 즉 $b(A) = \{b, e\}$이다.

Example 5.3: 유리수의 집합 $\mathbb{Q}$를 생각하자. $\mathbb{R}$의 모든 열린 부분집합은 유리수와 무리수를 모두 포함하므로, $\mathbb{Q}$의 내점이나 외점은 존재하지 않는다. 따라서 $\text{int}(\mathbb{Q}) = \varnothing$이고 $\text{int}(\mathbb{Q}^c) = \varnothing$이다. 따라서 $\mathbb{Q}$의 경계는 실수 전체의 집합, 즉 $b(\mathbb{Q}) = \mathbb{R}$이다.

위상 공간 $X$의 부분집합 $A$가 $X$에서 어디에서도 조밀하지 않다란, $A$의 폐포의 내부가 공집합인 것, 즉 $\text{int}(\bar{A}) = \varnothing$인 것이다.

Example 5.4: $\mathbb{R}$의 부분집합 $A = \{1, \frac{1}{2}, \frac{1}{3}, \frac{1}{4}, ...\}$를 생각하자. 앞서 언급했듯이 $A$는 정확히 하나의 극한점 $0$을 가진다. 따라서 $\bar{A} = \{0, 1, \frac{1}{2}, \frac{1}{3}, \frac{1}{4}, ...\}$이다. $\bar{A}$는 내점을 갖지 않음을 확인하라. 따라서 $A$는 $\mathbb{R}$에서 어디에서도 조밀하지 않다.

Example 5.5: $A$를 0과 1 사이의 유리점들로 이루어진 집합, 즉 $A = \{x : x \in \mathbb{Q}, 0 < x < 1\}$이라 하자. $A$의 내부는 공집합, 즉 $\text{int}(A) = \varnothing$임을 확인하라. 그러나 $A$는 $\mathbb{R}$에서 어디에서도 조밀하지 않은 것은 아니다. $A$의 폐포가 $[0, 1]$이므로,
$$\text{int}(\bar{A}) = \text{int}([0, 1]) = (0, 1)$$
은 공집합이 아니기 때문이다.