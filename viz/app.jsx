import React from "react";
import ReactDOM from "react-dom";
import Navbar from "react-bootstrap/Navbar";
import Nav from "react-bootstrap/Nav";
import Modal from "react-bootstrap/Modal";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import { json } from "d3";
import PDP from "alignment.js/prevent_default_patch";

import Superreads from "./superreads.jsx";

PDP(document);

import "./style.scss";


class Header extends React.Component {
  constructor(props) {
    super(props);
    this.state = { show: false };
  } 
  render() {
    return (<div>
      <Navbar bg="light">
        <Navbar.Brand>SuperSEAL</Navbar.Brand>
        <Nav.Item className="ml-auto">
          <Nav.Link onClick={() => this.setState({show: true})}>
            About
          </Nav.Link>
        </Nav.Item>
      </Navbar>
      <Container>
        {this.props.children}
      </Container>
      <Modal
        show={this.state.show}
        onHide={() => this.setState({show: false})}
      >
        <Modal.Header closeButton>
        </Modal.Header>
        <Modal.Body>
          <div style={{
              display: "flex",
              alignItems: "center",
              flexDirection: "column"
            }}
          >
            <img style={{width: "50%", height: "50%" }} src="/logo.png" />
            <h2>SuperSEAL</h2>
            <b>A viral quasispecies reconstruction toolkit</b>
            <p>Written by Stephen D. Shank, Ph. D.</p>
            <p>An <a href="http://lab.hyphy.org">
              Acme Computational Molecular Evolution
            </a> product</p>
            <a href="https://github.com/stephenshank/superseal">
              Documentation and source code
            </a>
          </div>
        </Modal.Body>
      </Modal>
    </div>);
  }
}

class App extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      superreads: null,
      viewing: "superreads"
    }
  }
  componentDidMount() {
    json("/superreads.json")
      .then(json => {
        json.number_of_sites = json.map(superread => superread.cv_end)
          .reduce((acc, curr) => Math.max(acc, curr), 0);
        this.setState({superreads: json});
      });
  }
  render() {
    switch (this.state.viewing) {
      case "superreads":
        return(<Header>
          <Superreads json={this.state.superreads} />
        </Header>);
    }
  }
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement("div"))
);
