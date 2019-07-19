import React, { Component } from "react";

class ErrorCorrection extends Component {
  constructor(props) {
    super(props);
    this.state = {
      error: null
    }
  }
  render() {
    return (<div>
      <h1>Error Correction</h1>
    </div>);
  }
}

export default ErrorCorrection;
