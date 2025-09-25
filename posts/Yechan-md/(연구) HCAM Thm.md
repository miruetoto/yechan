---
title: (연구) Hcam Thm
author: 신록예찬
date: 08/05/2025
draft: false
---

## 첫번째 버전

### 정보 이론적 관점

H-CAM의 계층적 분해는 정보 이론의 상호 정보량 분해와 유사한 구조를 가짐.

\begin{defn}[상호 정보량]
입력 이미지 $\mathbf{X}$와 클래스 라벨 $Y$ 사이의 상호 정보량은 다음과 같이 정의됨:
$$I(\mathbf{X}; Y) = H(Y) - H(Y|\mathbf{X})$$
여기서 $H(Y)$는 $Y$의 엔트로피이고, $H(Y|\mathbf{X})$는 $\mathbf{X}$가 주어졌을 때 $Y$의 조건부 엔트로피임. 이는 이미지 $\mathbf{X}$가 클래스 예측에 제공하는 정보량을 나타냄.
\end{defn}

\begin{thm}[H-CAM의 정보 분해 정리]
H-CAM에서 생성되는 각 계층 $\mathbf{X}^{(k)}$ (k=1,2,…,K)에 대해 다음이 성립함:
$$I(\mathbf{X}; Y) \geq \sum_{k=1}^{K} I(\mathbf{X}^{(k)}; Y)$$
여기서 등호는 각 계층이 서로 독립이고 정보 손실이 없을 때 성립함.
\end{thm}

\begin{proof}
H-CAM의 분해 과정에서 $\mathbf{X} = \sum_{k=1}^{K} \mathbf{X}^{(k)} + \mathbf{R}$로 표현할 수 있음. 여기서 $\mathbf{R}$은 잔여 성분임.

정보 이론의 부등식에 의해:
$$I(\mathbf{X}; Y) = I(\sum_{k=1}^{K} \mathbf{X}^{(k)} + \mathbf{R}; Y)$$
상호정보량의 체인규칙에 의하여 아래가 성립한다: 
$$I(\mathbf{X}^{(1)}, \mathbf{X}^{(2)}, \ldots, \mathbf{X}^{(K)}; Y) \geq \sum_{k=1}^{K} I(\mathbf{X}^{(k)}; Y)$$
마스킹 과정에서 정보 손실이 발생할 수 있으므로:
$$I(\mathbf{X}; Y) \geq I(\mathbf{X}^{(1)}, \mathbf{X}^{(2)}, \ldots, \mathbf{X}^{(K)}; Y)$$

따라서 연쇄 부등식에 의해 원하는 결과를 얻음.
\end{proof}

\begin{thm}[계층별 정보량 감소 정리]
H-CAM에서 계층 인덱스 $i$가 증가할수록 다음이 성립함:
$$I(\mathbf{X}^{(i)}; Y) \geq I(\mathbf{X}^{(i+1)}; Y)$$
또한 이는 CAM 분산의 감소와 밀접한 관련이 있음:
$$\sigma^2(\mathbf{M}^{(i)}) \geq \sigma^2(\mathbf{M}^{(i+1)})$$
\end{thm}

\begin{proof}
H-CAM의 마스킹 과정에서 $\mathbf{X}^{(i+1)}$은 $\mathbf{X}^{(i)}$에서 가장 중요한 특징을 제거한 잔여 부분임.

정보 처리 부등식(Data Processing Inequality)에 의해, 정보를 제거하는 과정은 상호 정보량을 감소시키거나 유지할 뿐 증가시킬 수 없음:
$$I(\mathbf{X}^{(i)}; Y) \geq I(\mathbf{X}^{(i+1)}; Y)$$

CAM 분산과의 관계에서, $\mathbf{M}^{(i)}$는 $\mathbf{X}^{(i)}$의 클래스별 중요도를 나타내므로:
$$\sigma^2(\mathbf{M}^{(i)}) \propto I(\mathbf{X}^{(i)}; Y)$$

따라서 정보량이 감소하면 CAM의 분산도 감소함.
\end{proof}

이러한 정보 이론적 분석은 H-CAM이 원본 이미지의 분류 정보를 계층별로 체계적으로 분해하며, 각 계층에서 추출되는 정보량이 단조 감소함을 보여줌. 또한 CAM 분산의 감소가 정보량 감소와 이론적으로 연결되어 있음을 확인할 수 있음.​​​​​​​​​​​​​​​​

## 정보 이론적 관점

H-CAM의 계층적 분해는 정보 이론의 상호 정보량 분해와 유사한 구조를 가짐.

\begin{defn}[상호 정보량]
입력 이미지 $\mathbf{X}$와 클래스 라벨 $Y$ 사이의 상호 정보량은 다음과 같이 정의됨:
$$I(\mathbf{X}; Y) = H(Y) - H(Y|\mathbf{X})$$
여기서 $H(Y)$는 $Y$의 엔트로피이고, $H(Y|\mathbf{X})$는 $\mathbf{X}$가 주어졌을 때 $Y$의 조건부 엔트로피임. 이는 이미지 $\mathbf{X}$가 클래스 예측에 제공하는 정보량을 나타냄.
\end{defn}

\begin{thm}[Data Processing Inequality]
확률변수 $X$, $Y$, $Z$가 마르코프 체인 $X \to Y \to Z$를 형성하면:
$$I(X; Z) \leq I(X; Y)$$
즉, 정보 처리 과정에서 상호 정보량은 감소하거나 유지될 뿐 증가하지 않음.
\end{thm}

\begin{defn}[정규화 정보 손실]
최대값 정규화 연산 $\mathcal{N}: \mathbb{R}^{H \times W \times C} \to \mathbb{R}^{H \times W \times C}$에 대해, 정규화로 인한 정보 손실을 다음과 같이 정의함:
$$\Delta I_{\text{norm}} = I(\mathbf{X}; Y) - I(\mathcal{N}(\mathbf{X}); Y)$$
이는 동적 범위 압축으로 인해 손실되는 분류 정보량을 나타냄.
\end{defn}

\begin{thm}[H-CAM의 완화된 정보 분해 정리]
H-CAM에서 각 계층 $k$에 대해 새로운 모델 $f^{(k)}$를 훈련할 때, 다음이 성립함:
$$I(\mathbf{X}; Y) \geq \mathbb{E}\left[\sum_{k=1}^{K} I(\mathbf{X}^{(k)}; Y^{(k)})\right] + \sum_{k=1}^{K} \Delta I_{\text{norm}}^{(k)} + \mathcal{R}_{\text{overlap}}$$
여기서:
\begin{itemize}
\item $Y^{(k)}$는 모델 $f^{(k)}$의 예측 출력
\item $\Delta I_{\text{norm}}^{(k)}$는 $k$번째 계층의 정규화로 인한 정보 손실
\item $\mathcal{R}_{\text{overlap}} \geq 0$는 계층 간 정보 중복량
\item 기댓값은 훈련 데이터 분포에 대해 취함
\end{itemize}
\end{thm}

\begin{proof}
**Step 1**: 각 계층에서 모델을 재훈련하므로, $k$번째 모델 $f^{(k)}$는 마스킹된 데이터 $\mathcal{D}^{(k)}$에 적응됨. 이때 분류 정보량은:
$$I(\mathbf{X}^{(k)}; Y^{(k)}) = I(\mathbf{X}_{\text{res}}^{(k-1)} \odot \mathbf{W}^{(k)}; Y^{(k)})$$

**Step 2**: 최대값 정규화 적용 시 정보 손실이 발생함:
$$I(\mathbf{X}_{\text{norm}}^{(k)}; Y^{(k)}) = I(\mathbf{X}^{(k)}; Y^{(k)}) - \Delta I_{\text{norm}}^{(k)}$$

정규화로 인한 손실은 다음과 같이 근사할 수 있음:
$$\Delta I_{\text{norm}}^{(k)} \approx \log\left(\frac{\text{range}(\mathbf{X}^{(k)})}{\max(\mathbf{X}^{(k)})}\right) \cdot H(Y^{(k)})$$

**Step 3**: H-CAM의 마스킹은 공간적으로 분리되지만 의미적으로 중복될 수 있는 영역들을 생성함. 따라서:
$$I(\mathbf{X}^{(1)}, \mathbf{X}^{(2)}, \ldots, \mathbf{X}^{(K)}; Y) \leq \sum_{k=1}^{K} I(\mathbf{X}^{(k)}; Y^{(k)}) - \mathcal{R}_{\text{overlap}}$$

**Step 4**: Data Processing Inequality와 마스킹 과정의 정보 손실을 고려하면:
$$I(\mathbf{X}; Y) \geq I(\mathbf{X}^{(1)}, \mathbf{X}^{(2)}, \ldots, \mathbf{X}^{(K)}; Y) + \sum_{k=1}^{K} \Delta I_{\text{norm}}^{(k)}$$

**Step 5**: 모델의 불완전한 학습과 확률적 변동성을 고려하여 기댓값을 취하면:
$$I(\mathbf{X}; Y) \geq \mathbb{E}\left[\sum_{k=1}^{K} I(\mathbf{X}^{(k)}; Y^{(k)})\right] + \sum_{k=1}^{K} \Delta I_{\text{norm}}^{(k)} + \mathcal{R}_{\text{overlap}}$$
\end{proof}

\begin{cor}[실용적 함의]
위 정리는 다음을 시사함:
\begin{enumerate}
\item H-CAM의 각 계층은 원본의 부분 정보를 포함하며, 정보의 총합은 원본을 초과하지 않음
\item 정규화로 인한 정보 손실 $\Delta I_{\text{norm}}^{(k)}$는 피할 수 없지만, 적절한 $\theta$ 선택으로 최소화 가능
\item 계층 간 정보 중복 $\mathcal{R}_{\text{overlap}}$는 H-CAM이 완전히 독립적인 특징을 추출하지 못함을 의미하지만, 이는 실제로 해석의 강건성을 높일 수 있음
\end{enumerate}
\end{cor}

이 완화된 정리는 H-CAM이 원본 이미지의 분류 정보를 **손실 있는 계층적 분해(lossy hierarchical decomposition)**로 재구성함을 보여줌. 각 계층은 서로 다른 중요도의 특징들을 포착하며, 정규화와 중복으로 인한 정보 손실에도 불구하고 전체적으로는 원본의 주요 분류 근거들을 체계적으로 드러냄. 이는 H-CAM이 완벽한 정보 보존이 아닌 **해석 가능한 근사(interpretable approximation)**를 제공한다는 실용적 관점을 이론적으로 뒷받침함.

[[Yechan 인덱스]]
