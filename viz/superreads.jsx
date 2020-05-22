import React from "react";
import { AxisLeft, AxisTop } from "d3-react-axis";
import { scaleLinear, range, max } from "d3";

import Button from "react-bootstrap/Button";
import Modal from "react-bootstrap/Modal";
import ReactJson from "react-json-view";
import BaseAlignment from "alignment.js/components/BaseAlignment";
import SequenceAxis from "alignment.js/components/SequenceAxis";
import BaseSiteStackedBarChart from "alignment.js/components/BaseSiteStackedBarChart";
import BaseSequenceBarPlot from "alignment.js/components/BaseSequenceBarPlot";
import ScrollBroadcaster from "alignment.js/helpers/ScrollBroadcaster";
import {
  nucleotide_color, nucleotide_colors
} from "alignment.js/helpers/colors";
import computeLabelWidth from "alignment.js/helpers/computeLabelWidth";
import css_grid_format from "alignment.js/helpers/format";


const valid_characters = ["A", "C", "G", "T", "-"];

function filter_and_count(superreads, weight_filter, index_filter) {
  const filtered_superreads = superreads.filter(superread => {
    const { weight, index } = superread,
      right_weight = weight >= weight_filter || 0,
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
      if (valid_characters.includes(character)) {
        counts[character][site_index] += superread.weight;
        counts.total[site_index] += superread.weight;
      }
    }
  });
  const sequence_data = filtered_superreads.map(superread => {
    const head_gaps = '-'.repeat(superread.cv_start),
      tail_gaps = '-'.repeat(superreads.number_of_sites - superread.cv_end),
      seq = head_gaps + superread.vacs + tail_gaps;
    return {
      seq: seq,
      header: 'superread-' + (superread.index)
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
    }),
    weights: filtered_superreads.map(sr => sr.weight)
  };
}



function Help() {
  return (<div>
    <p>This visualization displays superreads restricted to covarying sites,
    along with their weights and nucleotide frequencies.</p>
    <p><b>Wheel</b> while hovering over the alignment to scroll through it.</p>
    <p><b>Filter</b> superreads by their weight by inputting to the weight filter.</p>
    <p><b>Click</b> a superread label to highlight differences among other superreads.
    Click the same label to unhighlight.</p>
    <p><b>View</b> underlying JSON by clicking the `View JSON` button.</p>
  </div>);
}


class Superreads extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      weight_filter: 0,
      show_json: false,
      show_help: false,
      index_filter: null,
      highlight: null
    }
  }
  render() {
    if (!this.props.json) return null;
    const { sequence_data, counts, weights } = filter_and_count(
        this.props.json, this.state.weight_filter, this.state.index_filter
      ),
      { width, bar_height, height, site_size, label_padding, bar_width } = this.props,
      label_width = computeLabelWidth(sequence_data, label_padding),
      full_pixel_width = sequence_data[0].seq.length * site_size,
      full_pixel_height = sequence_data.length * site_size,
      base_alignment_width = width - label_width,
      base_alignment_height = height - bar_height,
      alignment_width = Math.min(full_pixel_width, base_alignment_width),
      alignment_height = Math.min(full_pixel_height, height - bar_height),
      frequency_scale = scaleLinear().domain([0, 1]).range([bar_height-1, 0]),
      weight_scale = scaleLinear().domain([0, max(weights)]).range([0, bar_width-1]),
      container_style = {
        display: "grid",
        gridTemplateColumns: css_grid_format([label_width, alignment_width, bar_width]),
        gridTemplateRows: css_grid_format([bar_height, alignment_height])
      },
      max_weight = this.props.json.reduce(
        (acc, curr) => Math.max(acc, curr.weight), 0
      ),
      label_fills = Array(sequence_data.length).fill("black"),
      site_color = this.state.highlight != null ?
        (mol, site, header) => {
          const chosen_superread = sequence_data[this.state.highlight];
          if (header == chosen_superread.header) return nucleotide_colors[mol];
          if (mol == "-") return nucleotide_colors["-"];
          return mol == chosen_superread.seq[site - 1]
            ? "#EEE"
            : nucleotide_colors[mol];
        } : nucleotide_color,
      molecule = this.state.highlight != null ?
        (mol, site, header) => {
          const chosen_superread = sequence_data[this.state.highlight];
          if (header == chosen_superread.header) return mol;
          if (mol == "-") return "-";
          return mol == chosen_superread.seq[site - 1] ? "." : mol;
        } : mol => mol,
      scroll_broadcaster = new ScrollBroadcaster({
        width: full_pixel_width,
        height: full_pixel_height,
        x_pad: base_alignment_width,
        y_pad: base_alignment_height,
        bidirectional: [
          "alignmentjs-alignment",
          "alignmentjs-stacked-bar",
          "alignmentjs-labels-div",
          "alignmentjs-bar"
        ]
      });
    if(this.state.highlight != null) {
      label_fills[this.state.highlight] = "red";
    }
    return (<div>
      <div className="toolbar">
        <span className="toolbar-item">
          Weight filter:
        </span>
        <input
          type="number"
          min="0"
          value={this.state.weight_filter}
          onChange={e => this.setState({
            weight_filter: Math.min(e.target.value, max_weight) || null}
          )}
          style={{width: 50}}
        />
        <Button
          onClick={()=>this.setState({show_json: true})}
          className="float-right"
        >
          View JSON
        </Button>
        <Button
          onClick={()=>this.setState({show_help: true})}
          className="float-right"
        >
          Help
        </Button>

      </div>
      <div id="alignmentjs-main-div" style={container_style}>
        <svg width={label_width} height={bar_height}>
          <text
            x={label_width-40}
            y={bar_height/2}
            textAnchor="middle"
            alignmentBaseline="middle"
            transform={`rotate(-90 ${label_width-40} ${bar_height/2})`}
          >
            Frequency
          </text>
          <AxisLeft
            scale={frequency_scale}
            transform={`translate(${label_width-1}, 0)`}
            tickValues={range(.1, 1, .1)}
          />
        </svg>
        <BaseSiteStackedBarChart
          width={alignment_width}
          height={bar_height}
          data={counts}
          scroll_broadcaster={scroll_broadcaster}
        />
        <svg width={bar_width} height={bar_height}>
          <text
            x={bar_width/2}
            y={bar_height - 30}
            textAnchor="middle"
            alignmentBaseline="middle"
          >
            Weight
          </text>
          <AxisTop
            scale={weight_scale}
            transform={`translate(0, ${bar_height-1})`}
            tickValues={weight_scale.ticks(3).slice(1)}
          />
        </svg>
        <SequenceAxis
          width={label_width}
          height={alignment_height}
          sequence_data={sequence_data}
          site_size={site_size}
          scroll_broadcaster={scroll_broadcaster}
          fill={label_fills}
          onClick={(label, i) => {
            if(this.state.highlight == i) {
              this.setState({highlight: null})
            } else {
              this.setState({highlight: i})
            }
          }}
        />
        <BaseAlignment
          sequence_data={sequence_data}
          width={alignment_width}
          height={alignment_height}
          site_color={site_color}
          site_size={this.props.site_size}
          molecule={molecule}
          scroll_broadcaster={scroll_broadcaster}
        />
        <BaseSequenceBarPlot
          data={weights}
          width={bar_width}
          height={alignment_height}
          scroll_broadcaster={scroll_broadcaster}
          scale={weight_scale}
        />
      </div>
      <Modal
        size="lg"
        show={this.state.show_json}
        onHide={() => this.setState({show_json: false})}
      >
        <Modal.Header closeButton>
          <Modal.Title>Superread JSON</Modal.Title>
        </Modal.Header>

        <Modal.Body>
          <div style={{height: 500, overflowY: "scroll"}}>
            <ReactJson
              src={this.props.json}
              displayDataTypes={false}
              collapsed={1}
            />
          </div>
        </Modal.Body>
      </Modal>
      <Modal
        show={this.state.show_help}
        onHide={() => this.setState({show_help: false})}
      >
        <Modal.Header closeButton>
          <Modal.Title>Help</Modal.Title>
        </Modal.Header>

        <Modal.Body>
          <Help />
        </Modal.Body>
      </Modal>
    </div>);
  }
}


Superreads.defaultProps = {
  label_padding: 10,
  left_bar_padding: 10,
  right_bar_padding: 20,
  site_size: 20,
  bar_height: 100,
  bar_width: 200,
  width: 960,
  height: 500,
  sender: "main",
};


export default Superreads;
