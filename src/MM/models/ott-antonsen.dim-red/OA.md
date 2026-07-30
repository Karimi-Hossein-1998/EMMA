# Ott-Antonsen Ansatz: Dimensionality Reduction

In this section, we will explore the Ott-Antonsen (OA) Ansatz and its application in reducing the dimensionality of systems described by partial differential equations (PDEs). The OA Ansatz is particularly useful in the context of kinetic theory and plasma physics, where it can help simplify the description of complex many-body systems.

The key idea behind the OA Ansatz is to reduce the infinite-dimensional distribution function of the system to a finite set of moments. This is achieved by assuming a specific functional form for the distribution function, which is assumed to be a symmetrical Fourier expansion. By doing so, and assuming a Lorenzian form for the distribution function of the intrinsic frequencies, one can further simplify the analysis and obtain explicit solutions for the reduced system. The OA ansatz also imposes another condition on the distribution function of the oscillators, which is essentially assuming that the higher moments can be obtained from the first moment.

The continuum limit applies on $ f(\theta, \omega; t) $:

$$
\displaystyle \frac{\partial }{\partial t} f + \frac{\partial}{\partial \theta} (f \dot{\theta}) = 0.
$$

Here $ f(\theta, \omega; t) $ is the distribution function of the system, which depends on the phase space variables $ \theta $ --phases--and $ \omega $--intrinsic frequencies.

## Kuramoto Model and Continuum Limit

The Kuramoto model is normally written as,

$$
\displaystyle \frac{\mathrm{d}}{\mathrm{d} t} \theta_j(t) = \omega_j + \frac{K}{N} \sum_{k=1}^{N} \sin\left(\theta_k\left(t\right)-\theta_j\left(t\right)\right)
$$

The order parameter is defined as,

$$
\displaystyle r\left(t\right) = \rho\left(t\right) e^{i \varphi\left(t\right)} = \frac{1}{N} \sum_{j=1}^{N} e^{i \theta_j \left(t\right)}.
$$

Which is evidently the first moment of the Fourier expansion basis i.e. $ r(t) = \left\langle e^{i \Theta(t)} \right\rangle $.

## The Continuous System

If we imagine that there are infinitely many oscillators, we can transition from the discrete Kuramoto model to a continuous description. In this limit, the summation becomes an integral over the phase space, and we can define a density function $f(\theta, \omega; t)$ that describes the distribution of oscillators in the space of phases $\theta$ and intrinsic frequencies $\omega$.

Hence for the continuity equation, the Kuramoto model and the order parameter we respectively have:

$$
\displaystyle \frac{\partial}{\partial t} f + \frac{\partial}{\partial \theta} \left(f \dot{\theta}\right) = 0,
$$

,

$$
\displaystyle \frac{\mathrm{d} \theta(t)}{\mathrm{d} t} = \omega + K \cdot \int_{-\infty}^{\infty} \int_{-\pi}^{\pi} f(\theta', \omega'; t) \sin(\theta' - \theta(t)) \mathrm{d}\theta' \mathrm{d}\omega',
$$

and

$$
\displaystyle r(t) = \int_{-\infty}^{\infty}\int_{-\pi}^{\pi} f(\theta, \omega; t) e^{i \theta} \mathrm{d}\theta \mathrm{d}\omega.
$$

Which can be interpreted as replacing $\frac{1}{N} \sum_{k=1}^{N}\cdots$ with $\int_{-\infty}^{\infty} \int_{-\pi}^{\pi} \cdots \mathrm{d}\theta \mathrm{d}\omega$.

## Ott-Antonsen Ansatz

Consider $f(\theta, \omega; t) = \rho(\theta, \omega; t) e^{i \theta}$, where $\rho(\theta, \omega; t)$ is a real-valued function. The Ott-Antonsen Ansatz assumes that the distribution function can be expressed in this form, which allows us to focus on the dynamics of the phase and the amplitude separately.

Consider $f(\theta,\omega;t)$ as the general form of the continuous system of Kuramoto oscillators, furthermore consider $f$ to be a product of the distribution function of the intrinsic frequencies $g(\omega)$ and an arbitrary (hopefully!) analytic function $p(\theta,\omega;t)$ so that it is possible to write it as a Fourier series:

$$
f(\theta, \omega; t) = g(\omega) \sum_{n=-\infty}^{\infty} a_n(\omega; t) e^{i n \theta},
$$

where the coefficients $a_n(\omega; t)$ are determined by the dynamics of the system. We also know that the function is ought to be symmetric, i.e. $a_n(\omega; t) = a_{-n}^*(\omega; t)$.

The final component that we need to consider is the normalization condition for the distribution function, which ensures that the total density of oscillators is conserved:

$$
\int_{-\infty}^{\infty} \int_{-\pi}^{\pi} f(\theta, \omega; t) \mathrm{d}\theta \mathrm{d}\omega = 1.
$$

which yields $a_0 = \frac{1}{2\pi}$. And by extending to the complex plane and assuming $g$ to be Lorenzian we get:

$$
\displaystyle r_l = \int_{-\infty}^{\infty} \int_{-\pi}^{\pi} f(\theta, \omega; t) e^{i l \theta} \mathrm{d}\theta \mathrm{d}\omega = a_{l}^{\ast}(\tilde{\omega};t).
$$

\bftext{OR}

$$
\left\langle \left(e^{i \Theta}\right)^l \right\rangle = \int_{-\infty}^{\infty} \int_{-\pi}^{\pi} f(\theta, \omega; t) e^{i l \theta} \mathrm{d}\theta \mathrm{d}\omega = a_{l}^{\ast}(\tilde{\omega};t).
$$

the *constant* $\tilde{\omega}$ is the average intrinsic frequency of the oscillator distribution or in the complex plane it represents the pole of the distribution $g(\omega)$.
For instance, $r(t) = \left\langle e^{i \Theta} \right\rangle = a_{1}^{\ast}(\tilde{\omega};t)$ relates the order parameter and the distribution function (the arbitrary function $p(\theta,\omega;t)$).
Then we can find the partial differential equation(s) governing the dynamics of the coefficient functions $a_n$.

$$
\displaystyle \frac{\partial }{\partial t}a_n(\omega;t) + i n \omega a_n(\omega;t) + n\cdot\frac{K}{2} \left( a_{n+1}\left(\omega;t\right) r\left(t\right) - a_{n-1}(\omega;t) r^{\ast}\left(t\right) \right) = 0.
$$

assuming $a_n\left(\omega;t\right)=\left(a\left(\omega;t\right)\right)^n ,\; \forall n>0$ we conclude,

$$
\displaystyle \frac{\partial }{\partial t}a(\omega;t) + i \omega a(\omega;t) + \cdot\frac{K}{2} \left( \left(a\left(\omega;t\right)\right)^2 r\left(t\right) - r^{\ast}\left(t\right) \right) = 0.
$$

and then by using $r\left(t\right)=a^{\ast}\left(\tilde{\omega};t\right)$,

$$
\displaystyle \frac{\mathrm{d}}{\mathrm{d} t} r\left(t\right) - i \left(\tilde{\omega}\right)^{\ast} r\left(t\right) + \frac{K}{2}\left(\left(r\left(t\right)\right)^2 r^{\ast}\left(t\right) - r\left(t\right)\right) = 0.
$$

OR

$$
\displaystyle \frac{\mathrm{d}}{\mathrm{d} t} r\left(t\right) = i \left(\tilde{\omega}\right)^{\ast} r\left(t\right) - \frac{K}{2}\left(\left(r\left(t\right)\right)^2 r^{\ast}\left(t\right) - r\left(t\right)\right).
$$

there is only one (complex) ordinary differential equation governing the dynamics of the our originally infinite-dimensional system of oscillators.
The Ott-Antonsen reduction has effectively transformed the problem into a lower-dimensional one, capturing the essential dynamics while discarding irrelevant degrees of freedom.
We can further eplore the reduced system by assuming a specific form for the distribution of the intrinsic frequencies, in other word $g\left(\omega\right)$. Let's say that we have:

$$
\displaystyle g\left(\omega\right) = \frac{\gamma}{\pi}\cdot\frac{1}{\left(\omega - \mu\right)^{2} + \gamma^{2}}.
$$

here $\tilde{\omega} = \mu\pm i \gamma$, but for if we confine $\omega$ into the $\R$ we have $\tilde{\omega} = \mu$. To solve the integral that relates $r\left(t\right)$ and $a\left(\omega;t\right)$ we need the former one.

Hence, we can express $r\left(t\right)$ in terms of $a\left(\omega;t\right)$ as follows:

$$
\displaystyle r\left(t\right) = a^{\ast}\left(\mu - i\gamma;t\right)
$$

**Note**: the integral is calculated on the lower half of the complex plane due to the intention of bounding the magnitude of $a$ or $r$ in time.

So we have:

$$
\displaystyle \frac{\mathrm{d}}{\mathrm{d} t} r\left(t\right) = \frac{\mathrm{d}}{\mathrm{d} t} a^{\ast}\left(\mu - i\gamma;t\right),
$$

then,

$$
\displaystyle \frac{\mathrm{d}}{\mathrm{d} t} r\left(t\right) = \left(- \gamma + i \mu \right) r\left(t\right) - \frac{K}{2}\left(\left(r\left(t\right)\right)^2 r^{\ast}\left(t\right) - r\left(t\right)\right),
$$

we then get the pair of equations that describe the dynamics of the magnitude (order parameter) $\rho\left(t\right)$ and the mean phase of the system $\varphi\left(t\right)$:

$$
\displaystyle \frac{\mathrm{d}}{\mathrm{d} t}\rho\left(t\right) = - \gamma \rho\left(t\right) - \frac{K}{2}\cdot \left(\left(\rho\left(t\right)\right)^3-\rho\left(t\right)\right).\\[1.0cm]
\frac{\mathrm{d}}{\mathrm{d} t}\varphi\left(t\right) = \mu.
$$
