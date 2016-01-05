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
#ifndef RX_DEMOD_P25_H
#define RX_DEMOD_P25_H

#include <gnuradio/analog/feedforward_agc_cc.h>
#include <gnuradio/digital/mpsk_receiver_cc.h>
#include <gnuradio/digital/diff_phasor_cc.h>
#include <gnuradio/blocks/complex_to_arg.h>
#include <gnuradio/blocks/multiply_const_ff.h>

#if 0
#include <gnuradio/analog/quadrature_demod_cf.h>
#include <gnuradio/filter/fir_filter.h>
#include <op25/fsk4_demod_ff.h>

#include <gnuradio/msg_queue.h>
#endif

#include <op25/fsk4_slicer_fb.h>
#include <op25/decoder_bf.h>
#include "dsp/resampler_xx.h"

#include <gnuradio/hier_block2.h>
#include <vector>


class rx_demod_p25;

typedef boost::shared_ptr<rx_demod_p25> rx_demod_p25_sptr;

/*! \brief Return a shared_ptr to a new instance of rx_demod_p25.
 *
 * This is effectively the public constructor.
 */
rx_demod_p25_sptr make_rx_demod_p25(float quad_rate, float audio_rate, float max_dev);


/*! \brief P25 demodulator.
 *  \ingroup DSP
 *
 * This class implements a P25 demodulator using OP25.
 *
 */
class rx_demod_p25 : public gr::hier_block2
{

public:
    rx_demod_p25(float quad_rate=48000.0, float audio_rate=48000.0, float max_dev=600.0); // FIXME: could be private
    ~rx_demod_p25();

private:
    /* GR blocks */
    /* CQPSK */
    gr::analog::feedforward_agc_cc::sptr d_agc;           /*! AGC/prescaler to for Costas loop. */
    gr::digital::mpsk_receiver_cc::sptr  d_clock;         /*! Clock/phase recovery. */
    gr::digital::diff_phasor_cc::sptr    d_diffdec;       /*! Differential decoder. */
    gr::blocks::complex_to_arg::sptr     d_to_float;      /*! Complex to float. */
    gr::blocks::multiply_const_ff::sptr  d_rescale;       /*! Rescaler to {-3, -1, 1, 3} levels. */
#if 0
    /* C4FM */
    gr::analog::quadrature_demod_cf::sptr d_quad;          /*! The quadrature demodulator block. */
    gr::fir_filter_fff::sptr              d_symbol_filter; /*! Symbol filter. */
    gr::op25::fsk4_demod_ff::sptr         d_demod_fsk4;    /*! C4FM/FSK4 demodulator block. */
#endif
    /* both */
    gr::op25::fsk4_slicer_fb::sptr      d_slicer;        /*! FSK4 slicer block. */
    gr::op25::decoder_bf::sptr          d_op25;          /*! P25 decoder block. */
    resampler_ff_sptr                   d_audio_resamp;  /*! Audio resampler. */

    /* other parameters */
    float  d_quad_rate;        /*! Quadrature rate. */
    float  d_audio_rate;       /*! Audio rate. */
    float  d_max_dev;          /*! Max deviation. */
#if 0
    gr::msg_queue::sptr d_queue; /*! Message queue to FSK4 demod. */
#endif

    /* FIR filter taps */
    std::vector<float> d_symbol_coeffs; /*! Symbol coefficients. */

};


#endif // RX_DEMOD_P25_H
