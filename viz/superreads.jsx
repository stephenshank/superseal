import React from "react";
import { AxisLeft } from "d3-react-axis";
import { scaleLinear, max, range } from "d3";

import BaseAlignment from "alignment.js/components/BaseAlignment";
import SiteAxis from "alignment.js/components/SiteAxis";
import SequenceAxis from "alignment.js/components/SequenceAxis";
import Placeholder from "alignment.js/components/Placeholder";
import BaseSiteStackedBarChart from "alignment.js/components/BaseSiteStackedBarChart";
import ScrollBroadcaster from "alignment.js/helpers/ScrollBroadcaster";
import { nucleotide_color, nucleotide_text_color } from "alignment.js/helpers/colors";
import computeLabelWidth from "alignment.js/helpers/computeLabelWidth";
import css_grid_format from "alignment.js/helpers/format";


function filter_and_count(superreads, weight_filter, index_filter) {
  const filtered_superreads = superreads.filter(superread => {
    const { weight, index } = superread,
      right_weight = weight_filter ?
        weight >= weight_filter[0] && weight <= weight_filter[1] :
        true,
      right_index = index_filter ?
        index >= index_filter[0] && index <= index_filter[1] :
        true;
      return right_weight && right_index;
    });
    const counts = {
        A: range(superreads.number_of_sites).fill(0),
        C: range(superreads.number_of_sites).fill(0),
        G: range(superreads.number_of_sites).fill(0),
        T: range(superreads.number_of_sites).fill(0),
        total: range(superreads.number_of_sites).fill(0)
      };
  filtered_superreads.forEach(superread => {
    for(let i=0; i < superread.vacs.length; i++) {
      let character = superread.vacs[i],
        site_index = superread.cv_start + i;
      counts[character][site_index] += superread.weight;
      counts.total[site_index] += superread.weight;
    }
  });
  const sequence_data = superreads.map((superread, i) => {
    const head_gaps = '-'.repeat(superread.cv_start),
      tail_gaps = '-'.repeat(superreads.number_of_sites - superread.cv_end),
      seq = head_gaps + superread.vacs + tail_gaps;
    return {
      seq: seq,
      header: 'superread-' + (i+1)
    };
  });
  sequence_data.number_of_sequences = sequence_data.length;
  sequence_data.number_of_sites = sequence_data[0].length;
  return {
    sequence_data: sequence_data,
    counts: range(superreads.number_of_sites).map((x, i) => {
      return [
        counts.A[i] / counts.total[i],
        counts.C[i] / counts.total[i],
        counts.G[i] / counts.total[i],
        counts.T[i] / counts.total[i],
      ];
    })
  };
}

class Superreads extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      weight_filter: null,
      index_filter: null
    }
  }
  render() {
    if (!this.props.json) return null;
    const { sequence_data, counts } = filter_and_count(
        this.props.json, this.state.weight_filter, this.state.index_filter
      ),
      { width, bar_height, height, site_size, label_padding } = this.props,
      label_width = computeLabelWidth(sequence_data, label_padding),
      full_pixel_width = sequence_data[0].seq.length * site_size,
      full_pixel_height = sequence_data.length * site_size,
      base_alignment_width = width - label_width,
      base_alignment_height = height - bar_height,
      alignment_width = Math.min(full_pixel_width, base_alignment_width),
      alignment_height = Math.min(full_pixel_height, height - bar_height),
      container_style = {
        display: "grid",
        gridTemplateColumns: css_grid_format([label_width, alignment_width]),
        gridTemplateRows: css_grid_format([bar_height, alignment_height])
      },
      scroll_broadcaster = new ScrollBroadcaster({
        width: full_pixel_width,
        height: full_pixel_height,
        x_pad: base_alignment_width,
        y_pad: base_alignment_height,
        bidirectional: [
          "alignmentjs-alignment",
          "alignmentjs-stacked-bar",
          "alignmentjs-labels-div"
        ]
      });

    return (<div id="alignmentjs-main-div" style={container_style}>
      <Placeholder width={label_width} height={bar_height} />
      <BaseSiteStackedBarChart
        width={alignment_width}
        height={bar_height}
        data={counts}
        scroll_broadcaster={scroll_broadcaster}
      />
      <SequenceAxis
        width={label_width}
        height={alignment_height}
        sequence_data={sequence_data}
        site_size={site_size}
        scroll_broadcaster={scroll_broadcaster}
      />
      <BaseAlignment
        sequence_data={sequence_data}
        width={alignment_width}
        height={alignment_height}
        site_color={this.props.site_color}
        text_color={this.props.text_color}
        site_size={this.props.site_size}
        molecule={this.props.molecule}
        scroll_broadcaster={scroll_broadcaster}
      />
    </div>);
  }
}


Superreads.defaultProps = {
  site_color: nucleotide_color,
  text_color: nucleotide_text_color,
  label_padding: 10,
  left_bar_padding: 10,
  right_bar_padding: 20,
  site_size: 20,
  bar_height: 100,
  width: 960,
  height: 500,
  sender: "main",
  molecule: mol => mol
};


export default Superreads;
