# The Modular/Hierarchical Structure

In this section, we will discuss the modular and hierarchical structure of the Ott-Antonsen Ansatz and its implications for the dynamics of the system. The modular structure refers to the way in which the system can be decomposed into smaller, interacting modules, each of which can be analyzed independently. This is particularly useful in the context of complex systems, where the interactions between different components can give rise to emergent behavior.

## Kuramoto model (for connected communities)

In the context of the Kuramoto model, we can think of each oscillator as belonging to a specific community or module. These communities can interact with each other, leading to complex dynamics that are not present in a simple mean-field description. The modular structure allows us to analyze the behavior of each community independently, while still accounting for the interactions between them.

The formulation of it is written as,

$$
\displaystyle \frac{\mathrm{d}}{\mathrm{d} t} \theta_{j}^{\left(c\right)}\left(t\right) = \omega_{j} + \sum_{c' = 1}^{\left(c\right)} \frac{K_{c,c'}}{N_{c'}} \frac{N_{c'}}{N} \sum_{k \in c'} \sin\left(\theta_{k}^{c'}\left(t\right) - \theta_{j}^{\left(c\right)}\left(t\right)\right)
$$

here the indices $c$ and $c'$ are used to denote different communities or modules within the system. The index $c$ refers to the community to which the oscillator $j$ belongs, while the index $c'$ refers to a different community that is interacting with community $c$. This notation allows us to keep track of the interactions between different modules in the system. Also it is worth mentioning that $N = \sum_{c} N_{c}$, where $N_{c}$ is the number of oscillators in community $c$, also $\frac{K_{c,c'}}{N_{c'}} \frac{N_{c'}}{N}$ can be written as $\frac{K_{c,c'}}{N_{c'}} \eta_{c'}$ and refactor $K_{c,c'}\eta_{c'}$ so that it is equal to the new set of coupling strength. The $\eta_{c}$ is the fraction of the oscillators in community $c$.

Considering the explanations, the definition of the order parameter will need some adjustment. We need to consider the order parameter for each individual community and each compartment in the next layer (if the structure is hierarchical with hierarchy level greater than $1$). So, we would have,

$$
\displaystyle r^{\left(c\right)}\left(t\right) = \rho^{\left(c\right)}\left(t\right) e^{i \varphi^{\left(c\right)}\left(t\right)} = \frac{1}{N_c} \sum_{j=1}^{N_c} e^{i\theta_{j}^{\left(c\right)}\left(t\right)}.
$$

which is also,

$$
\displaystyle r^{\left(c\right)}\left(t\right) = \left\langle e^{i\Theta^{\left(c\right)}\left(t\right)} \right\rangle. \; \forall c \in \mathcal{C},
$$

this can be for the lowest level of hierarchy (or in another word, finest separation of communities) when $\mathcal{C} \equiv \{1,2,\cdots,C\}$, or can be representative of the whole system when $\mathcal{C} \equiv \{1\}$ and anything in between.

It should be clear that the global order parameter can be written in terms of the order parameters of the lower layers, $r\left(t\right)= \sum_{c \in \mathcal{C}} \eta_{c} r^{\left(c\right)}\left(t\right)$.

Given that we already dealt with one community and how we reduce it into a low-dimensional system we shall approach the dimensionality reduction for this system next.

## OA Ansatz For Communities

In order to apply the Ott-Antonsen Ansatz to each community, we can assume a similar form for the distribution function within each community:

$$
\displaystyle f^{\left(c\right)}(\theta, \omega; t) = \rho^{\left(c\right)}(\theta, \omega; t) e^{i \theta},
$$

and the rest follows the same pattern as before, leading to a set of coupled equations for the order parameters of each community.

$$
\displaystyle \frac{\mathrm{d}}{\mathrm{d} t} r^{\left(c\right)}\left(t\right) = \left( -\gamma^{\left(c\right)} + i \mu^{\left(c\right)}\right) r^{\left(c\right)}\left(t\right) - {\Large\sum_{c' = 1}^{C}} \left\{ K_{c,c'} \left(\left( r^{\left(c\right)}\left(t\right) \right)^{2} \left(r^{\left(c'\right)}\left(t\right)\right)^{\ast} - r^{\left(c'\right)}\left(t\right) \right) \right\}.
$$

We can further expand the equations of dynamics of the reduced system into two couple system of equations for amplitudes (orders) and the phases (expect phase of each community of oscillators), which is:

$$
\displaystyle \frac{\mathrm{d}}{\mathrm{d} t} \rho^{\left(c\right)}\left(t\right) = -\gamma^{\left(c\right)} \rho^{\left(c\right)}\left(t\right) - {\Large\sum_{c' = 1}^{C}} \left\{ K_{c,c'} \left(\left( \rho^{\left(c\right)}\left(t\right) \right)^{2} - 1 \right)\rho^{\left(c'\right)}\left(t\right)\cos\left(\varphi^{\left(c\right)}\left(t\right) - \varphi^{\left(c'\right)}\left(t\right)\right) \right\},\\[1.0cm]
\frac{\mathrm{d}}{\mathrm{d} t} \varphi^{\left(c\right)}\left(t\right) = \mu^{\left(c\right)} - {\Large\sum_{c' = 1}^{C}} \left\{ K_{c,c'} \left(\left( \rho^{\left(c\right)}\left(t\right) \right)^{2} + 1 \right)\frac{\rho^{\left(c'\right)}\left(t\right)}{\rho^{\left(c\right)}\left(t\right)}\sin\left(\varphi^{\left(c\right)}\left(t\right) - \varphi^{\left(c'\right)}\left(t\right)\right) \right\}.
$$

And this completely describes the dynamics of the order parameters for each community in the system.
