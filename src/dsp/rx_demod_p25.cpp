/* -*- c++ -*- */
/*
 * Copyright 2012-2016 Joshua Roys KK4AFZ.
 *
 * Gqrx is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Gqrx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Gqrx; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */
#include <gnuradio/io_signature.h>
#include <dsp/rx_demod_p25.h>
#include <algorithm>

#define SYMBOL_RATE 4800.0

/* Create a new instance of rx_demod_p25 and return a boost shared_ptr. */
rx_demod_p25_sptr make_rx_demod_p25(float quad_rate, float audio_rate, float max_dev)
{
    return gnuradio::get_initial_sptr(new rx_demod_p25(quad_rate, audio_rate, max_dev));
}

static const int MIN_IN = 1;  /* Mininum number of input streams. */
static const int MAX_IN = 1;  /* Maximum number of input streams. */
static const int MIN_OUT = 1; /* Minimum number of output streams. */
static const int MAX_OUT = 1; /* Maximum number of output streams. */

rx_demod_p25::rx_demod_p25(float quad_rate, float audio_rate, float max_dev)
    : gr::hier_block2 ("rx_demod_p25",
                      gr::io_signature::make (MIN_IN, MAX_IN, sizeof (gr_complex)),
                      gr::io_signature::make (MIN_OUT, MAX_OUT, sizeof (float))),
    d_quad_rate(quad_rate),
    d_audio_rate(audio_rate),
    d_max_dev(max_dev)
{
    float theta, loop_bw, fmax, fmin, mu, gain_mu, omega, gain_omega, omega_rel;
    float sl_arr[4] = {-2.0, 0.0, 2.0, 4.0};
    std::vector<float> slice_levels(sl_arr, sl_arr + sizeof(sl_arr) / sizeof(float));

    /* AGC/prescale */
    d_agc = gr::analog::feedforward_agc_cc::make(16, 1.0);

    /* Clock/phase recovery */
    theta = M_PI / 4.0; /* PI/4 DQPSK */
    loop_bw = 2.0 * M_PI / 100.0;
    fmax = 0.025;
    fmin = -fmax;
    mu = 0.0;
    gain_mu = 0.05;
    omega = d_quad_rate / SYMBOL_RATE;
    gain_omega = omega * omega / 4;
    omega_rel = 0.005; /* allow 0.5% variation in omega */
    d_clock = gr::digital::mpsk_receiver_cc::make(4, theta, loop_bw, fmin, fmax, mu, gain_mu, omega, gain_omega, omega_rel);

    /* Differential decoding */
    d_diffdec = gr::digital::diff_phasor_cc::make();

    /* Complex to float */
    d_to_float = gr::blocks::complex_to_arg::make();

    /* Rescale to {-3, -1, 1, 3} */
    d_rescale = gr::blocks::multiply_const_ff::make(1.0 / (M_PI / 4.0));

    /* FSK4 slicer */
    d_slicer = gr::op25::fsk4_slicer_fb::make(slice_levels);

    /* P25 decoder */
    d_op25 = gr::op25::decoder_bf::make();

    /* Audio resampler */
    d_audio_resamp = make_resampler_ff(d_audio_rate / 8000.0);

    connect(self(), 0, d_agc, 0);
    connect(d_agc, 0, d_clock, 0);
    connect(d_clock, 0, d_diffdec, 0);
    connect(d_diffdec, 0, d_to_float, 0);
    connect(d_to_float, 0, d_rescale, 0);
    connect(d_rescale, 0, d_slicer, 0);
    connect(d_slicer, 0, d_op25, 0);
    connect(d_op25, 0, d_audio_resamp, 0);
    connect(d_audio_resamp, 0, self(), 0);

#if 0
    float gain;
    float sl_arr[4] = {-2.0, 0.0, 2.0, 4.0};
    int samples_per_symbol;
    std::vector<float> slice_levels(sl_arr, sl_arr + sizeof(sl_arr) / sizeof(float));

    /* demodulator gain */
    gain = d_quad_rate / (2.0 * M_PI * d_max_dev);

    //std::cout << "G: " << gain << std::endl;

    /* demodulator */
    d_quad = gr_make_quadrature_demod_cf(gain);

    /* Symbol filter */
    samples_per_symbol = int(d_quad_rate / SYMBOL_RATE);
    d_symbol_coeffs.resize(samples_per_symbol);
    fill(d_symbol_coeffs.begin(), d_symbol_coeffs.end(), (1.0 / samples_per_symbol));
    d_symbol_filter = gr_make_fir_filter_fff(1, d_symbol_coeffs);

    /* C4FM/FSK4 demodulator */
    d_queue = gr_make_msg_queue();
    d_demod_fsk4 = op25_make_fsk4_demod_ff(d_queue, d_quad_rate, SYMBOL_RATE);

    /* FSK4 slicer */
    d_slicer = op25_make_fsk4_slicer_fb(slice_levels);

    /* P25 decoder */
    d_op25 = op25_make_decoder_bf();

    /* Audio resampler */
    d_audio_resamp = make_resampler_ff(d_audio_rate / 8000.0);

    /* connect block */
    connect(self(), 0, d_quad, 0);
    connect(d_quad, 0, d_symbol_filter, 0);
    connect(d_symbol_filter, 0, d_demod_fsk4, 0);
    connect(d_demod_fsk4, 0, d_slicer, 0);
    connect(d_slicer, 0, d_op25, 0);
    connect(d_op25, 0, d_audio_resamp, 0);
    connect(d_audio_resamp, 0, self(), 0);
#endif

}


rx_demod_p25::~rx_demod_p25 ()
{

}
