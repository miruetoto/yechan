---
title: (공부) Big-O, Little-o
author: 신록예찬
date: 08/14/2025
draft: false
---

| 번호  | 문제                                                 | 정답  |
| --- | -------------------------------------------------- | --- |
| 1   | $5n^2 + 3n + 1 = O(n^2)$                           | T   |
| 2   | $n = O(n^2)$                                       | T   |
| 3   | $n^2 = O(n)$                                       | F   |
| 4   | $\log n = O(n)$                                    | T   |
| 5   | $n^{0.5} = O(n^{0.6})$                             | T   |
| 6   | $n^{0.6} = O(n^{0.5})$                             | F   |
| 7   | $1 = O(\log n)$                                    | T   |
| 8   | $\log n = O(1)$                                    | F   |
| 9   | $n\log n = O(n^2)$                                 | T   |
| 10  | $n^2 = O(n^2)$                                     | T   |
| 11  | $n = o(n^2)$                                       | T   |
| 12  | $n^2 = o(n)$                                       | F   |
| 13  | $\log n = o(n)$                                    | T   |
| 14  | $\sqrt{n} = o(n)$                                  | T   |
| 15  | $n^{0.99} = o(n)$                                  | T   |
| 16  | $n^{1.01} = o(n)$                                  | F   |
| 17  | $1 = o(\log n)$                                    | T   |
| 18  | $\log n = o(1)$                                    | F   |
| 19  | $n = o(n)$                                         | F   |
| 20  | $n = O(n)$                                         | T   |
| 21  | $n = O_p(n)$                                       | T   |
| 22  | $n = o_p(n^2)$                                     | T   |
| 23  | $n^2 = o_p(n)$                                     | F   |
| 24  | $\sqrt{n} = O_p(1)$                                | F   |
| 25  | $n^{-1/2} = o_p(1)$                                | T   |
| 26  | $n^2 + 3n = \Theta(n^2)$                           | T   |
| 27  | $3n^2 = \Theta(n^2)$                               | T   |
| 28  | $n^2 = \Theta(n)$                                  | F   |
| 29  | $\log n = \Theta(\log_{10} n)$                     | T   |
| 30  | $n\log n = \Theta(n)$                              | F   |
| 31  | $n^{0.5} = \Theta(\sqrt{n})$                       | T   |
| 32  | $5n + 3\log n = \Theta(n)$                         | T   |
| 33  | $n^2 = \Theta(n^2)$                                | T   |
| 34  | $n = \Theta(n^2)$                                  | F   |
| 35  | $n^2 + n = \Theta(n)$                              | F   |
| 36  | $\log n = \Theta(\sqrt{\log n})$                   | F   |
| 37  | $n^{0.5} = \Theta(n^{0.6})$                        | F   |
| 38  | $n\log n = \Theta(n\log_{10} n)$                   | T   |
| 39  | $n\log n = \Theta(n\log^2 n)$                      | F   |
| 40  | $\sqrt{n} + \log n = \Theta(\sqrt{n})$             | T   |
| 41  | $n^3 + n^2 + 5 = \Theta(n^3)$                      | T   |
| 42  | $n^{1.5} + n = \Theta(n^2)$                        | F   |
| 43  | $n^2 = \Theta(n^2 + n)$                            | T   |
| 44  | $n^2 + n = \Theta(n^2 + \sqrt{n})$                 | T   |
| 45  | $\sqrt{n} = \Theta(n^{0.5})$                       | T   |
| 46  | $n^{0.5} = \Theta(n^{0.5})$                        | T   |
| 47  | $n^2\log n = \Theta(n^2)$                          | F   |
| 48  | $n^2\log n = \Theta(n^2\log_{10} n)$               | T   |
| 49  | $n^2 + 3n + 1 = \Theta(n^3)$                       | F   |
| 50  | $n^2 + n\log n = \Theta(n^2)$                      | T   |
| 51  | $n^2 + n \sim n^2$                                 | T   |
| 52  | $3n^2 \sim n^2$                                    | F   |
| 53  | $n + \log n \sim n$                                | T   |
| 54  | $n^3 + 5n^2 \sim n^3$                              | T   |
| 55  | $n^3 + n \sim n^2$                                 | F   |
| 56  | $\log n + 1 \sim \log n$                           | T   |
| 57  | $\sqrt{n} + n \sim n$                              | T   |
| 58  | $n(1 + 1/n) \sim n$                                | T   |
| 59  | $n(1 + 1/n)^2 \sim n$                              | T   |
| 60  | $n^2 + o(n^2) \sim n^2$                            | T   |
| 61  | $n^2(1 + o(1)) \sim n^2$                           | T   |
| 62  | $n^2 + \Theta(n) \sim n^2$                         | F   |
| 63  | $n^2 + O(n) \sim n^2$                              | F   |
| 64  | $n^2 + 2n^2 \sim n^2$                              | F   |
| 65  | $n + n^{0.5} \sim n$                               | T   |
| 66  | $n + 3n \sim n$                                    | F   |
| 67  | $n\log n + n \sim n\log n$                         | T   |
| 68  | $n\log n + \log n \sim n\log n$                    | T   |
| 69  | $n\log n + n^2 \sim n^2$                           | T   |
| 70  | $n^2 + n\log n \sim n^2$                           | T   |
| 71  | $n^2 + \log n \sim n^2$                            | T   |
| 72  | $\log(n^2) \sim 2\log n$                           | T   |
| 73  | $\log_{10} n \sim \log n$                          | T   |
| 74  | $\sqrt{n^2 + n} \sim n$                            | T   |
| 75  | $\sqrt{n^2 + n} \sim n + 0.5$                      | F   |
| 76  | $n^3 = o(n^4)$                                     | T   |
| 77  | $n^3 = O(n^4)$                                     | T   |
| 78  | $n^3 = \Theta(n^4)$                                | F   |
| 79  | $n^3 + 5n = \Theta(n^3)$                           | T   |
| 80  | $n^3 + 5n \sim n^3$                                | T   |
| 81  | $\log n = o(n^\epsilon) \ \forall \epsilon>0$      | T   |
| 82  | $n^\epsilon = o(\log n) \ \forall \epsilon>0$      | F   |
| 83  | $n\log n = o(n^{1+\epsilon}) \ \forall \epsilon>0$ | T   |
| 84  | $n^2\log n = O(n^3)$                               | T   |
| 85  | $n^2\log n = o(n^3)$                               | T   |
| 86  | $n^2\log n = \Theta(n^2)$                          | F   |
| 87  | $n^2\log n = \Theta(n^2\log_{10} n)$               | T   |
| 88  | $n^{1.5} = o(n^2)$                                 | T   |
| 89  | $n^{1.5} = O(n^2)$                                 | T   |
| 90  | $n^{1.5} = \Theta(n^2)$                            | F   |
| 91  | $n(1 + o(1)) = \Theta(n)$                          | T   |
| 92  | $n(1 + o(1)) \sim n$                               | T   |
| 93  | $n(1 + O(1)) = \Theta(n)$                          | T   |
| 94  | $n(1 + O(1)) \sim n$                               | F   |
| 95  | $n + o(n) = \Theta(n)$                             | T   |
| 96  | $n + o(n) \sim n$                                  | T   |
| 97  | $n + O(n) = \Theta(n)$                             | T   |
| 98  | $n + O(n) \sim n$                                  | F   |
| 99  | $n + \Theta(n) = \Theta(n)$                        | T   |
| 100 | $n + \Theta(n) \sim n$                             | F   |