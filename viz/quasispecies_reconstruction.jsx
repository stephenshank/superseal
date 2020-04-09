import React, { Component } from "react";
import SequenceBarChart from "alignment.js/SequenceBarChart";
import fastaParser from "alignment.js/helpers/fasta";
import axios from "axios";


function QuasispeciesAlignment(props) {
  if (!props.fasta) return null;
  const sequence_data = fastaParser(props.fasta),
    data = sequence_data.map(
      record => +record.header.split('_')[1].split('-')[1],
    );
  sequence_data.forEach(record => {
    record.header = record.header.split('_')[0]
  });
  return (<SequenceBarChart
    sequence_data={sequence_data}
    data={data}
    label="Quasispecies Frequency"
  />);
}

class QuasispeciesReconstruction extends Component {
  constructor(props) {
    super(props);
    this.state = {
      error: null,
      fasta: null
    }
  }
  componentDidMount() {
    axios.get('/api/quasispecies.fasta')
      .then(response => {
        this.setState({fasta: response.data})
      }).catch(error => {
        this.setState({error:true});
      });
  }
  render() {
    return (<div>
      <h1>Quasispecies Reconstruction</h1>
      <QuasispeciesAlignment fasta={this.state.fasta} />
    </div>);
  }
}

export default QuasispeciesReconstruction;
