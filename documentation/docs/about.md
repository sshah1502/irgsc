### About irgsctool
__irgsctool__ to generate the catalog of NIR guide stars for Adaptive Optics (AO) observations of the Thirty Meter Telescope (TMT). The module computes the NIR magnitudes of the optical sources in the PANSTARRS data.
### Motivation to generate irgsctool
<p style="text-align: justify;">The performance of any ground-based optical/near-infrared (NIR) telescope is affected by the turbulence in the atmosphere of Earth. When the light from a distant astronomical source passes through Earth's turbulent atmosphere, it distorts the wavefront of the light. These distortions make the science images appear fuzzy/blurry. To improve the performance of the ground-based optical/NIR telescopes by compensating for the effects of wavefront distortions, astronomers use a technique known as Adaptive Optics (AO). An AO system tries to correct the distortions using a WaveFront Sensor (WFS), which takes some of the astronomical light, a deformable mirror that lies in the optical path, and a computer that receives input from the detector (Refer Figure: (1)).</p>

<center>
<figure>
    <img src="/img/aos.png" width="600" height="500"
         alt="AOS">
    <figcaption><em><strong>Figure 1:</strong></em> A conceptual image of the working of an AO system.</figcaption>
</figure>
</center>

<p style="text-align: justify;">A WFS is a high-speed camera that detects the deformity in the incoming light several thousand times per second. After measuring the deformity, the WFS sends feedback to the system-controlling computer, which then changes the shape of the deformable mirror so that the light eventually becomes distortion-free. This distortion-free light is then fed to the science instruments in the telescopes to obtain high-quality images with a spatial resolution close to the theoretical diffraction limit. A science target to be studied is often too faint or extended to be used as a reference for measuring the shape of the incident wave-fronts. Instead, a nearby brighter guide star can be used for distortion correction. But, sufficiently bright stars are not available in all parts of the sky, which significantly limits the usefulness of natural guide star AO. This limitation can be overcome by creating an artificial guide star by shining a laser into the mesosphere and creating an asterism. However, much fainter natural reference stars are still required for image position or tip/tilt correction.</p>
<p style="text-align: justify;">The Thirty Meter Telescope (TMT) is one of the largest optical and near-infrared (NIR) ground-based telescopes to be built, and the first light is expected in the next decade. The first light NIR instruments of TMT will be assisted by a multi-conjugate AO instrument, known as Narrow Field Infrared Adaptive Optics System (NFIRAOS) (see Figure: (2)). NFIRAOS will have a field of view (FOV) of 2-arcmin diameter and will use laser guide stars (LGS) for distortion correction. However, it will require three or more natural guide stars for tip/tilt correction. A catalog of guide stars will thus be a critical resource for TMT operations. It will enable efficient planning and observations, fulfilling a role similar to that of the Guide Star Catalogs I and II, which were created to allow for the acquisition and control of the Hubble Space Telescope. The TMT Infrared Guide Star Catalog (TMT-IRGSC) should be a star catalog consisting of point sources with NIR (J, H, Ks) magnitudes as faint as 22 mag in the J band in the Vega system covering the entire TMT-observable sky. No catalog currently exists with objects as faint as J_Vega = 22 mag over the entire TMT observable sky to be used as a guide star catalog. Hence it is essential to develop this catalog by computing the expected NIR magnitudes of stellar sources identified in various optical sky surveys using their optical magnitudes.</p>

<justify>
<figure>
    <img src="/img/nfiraos.jpg" width="600" height="500"
         alt="nfiraos">
    <figcaption><em><strong>Figure 2:</strong></em> A rendered image of the facility AO system on TMT - NFIRAOS.</figcaption>
</figure>
</justify>